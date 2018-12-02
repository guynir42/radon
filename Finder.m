classdef Finder < handle
% Streak finding tool. 
% This object can be used in two ways:
%   (1) Give it some images, using obj.input(images, ...)
%       This will do all the internal calculations needed to find streaks
%       inside each frame. Optional arguments to "input" are:
%           -variance: give a scalar (average) variance *or* a variance map. 
%            At some point we may also allow giving a 3D matrix of
%            variance maps for each image in the batch. 
%           -psf: give the point spread function as a scalar width (in this 
%            case the PSF is a 2D Gaussian with sigma equal to the input 
%            scalar) *or* a map for the image. PSF can be 3D matrix with a 
%            number of slices as the number of images. 
%           -batch number: housekeeping parameter. Helps keep track of where 
%            each streak was found.
%           -filename: housekeeping parameter. Helps keep track of where 
%            each streak was found.
%
%   (2) Give the finder object as an optional argument to the "frt"
%       function, so it does its magic while the FRT is running. 
%       In this case you should make sure all parameters of the finder are
%       set correctly before starting the FRT on the image. 
%       You can construct a new finder in the function call itself, using
%       the constructor inputs to set the finder parameters, e.g.
%       [R,f] = radon.frt(im, 'finder', Finder('median', 1, 'short', 0,...));
%       In this case you must take 2 arguments from the "frt" function, the
%       second argument being the handle to the newly constructed finder. 
% 
% SWITCHES AND SETTINGS (see comments in properties block for more details)
%   -Image pre-processing: use_subtract_mean, use_conv, use_crop_image, crop_size.
%   -Search options: use_short, min_length, use_recursive, recursion_depth. 
%   -Post processing: use_exclude, exclude_dx, exclude_dy
%   -memory considerations: use_save_images, use_clear_memory.
%   -housekeeping: filename, batch_num, frame_num.
%   -display defaults: show_bit, display_index, display_which, line_offset,
%    rect_size, display_monochrome.
%
% NOTES: -the Radon image that is saved in each streak is already normalized 
%         by the variance map. 
%        -If you cut the stars out, replace them with NaNs, then subtract 
%         the mean of the image, then replace the NaNs with zeros.
%
% A Graphic User Interface (GUI) for this class is also included in the 
% sub-package +gui. It is invoked using obj.makeGUI. 

    properties(Transient=true) % GUI, variance maps
        
        gui@radon.gui.FinderGUI; % a GUI object for this class
        
    end
    
    properties % objects
        
        last_streak@radon.Streak; % the last (single) streak to be detected
        streaks@radon.Streak; % a vector of all streaks detected in the latest call to "input" 
        prev_streaks@radon.Streak; % a vector of all streaks since the last call to "reset" (all streaks per run)
        
        timing_data@util.time.TimingData; % a wrapper of a dictionary that tracks how much time is spent in each phase of the calculations
                
    end
        
    properties % inputs/outputs
                
        % images:
        input_images; % if we are given the image in real space
        im_size; % can be input as scalar, will be output as two-element vector.
        
        % input var/psf:
        input_var; % original variance map/scalar 
        input_psf; % original PSF as given (or width parameter "sigma" of PSF)        
        
        % outputs:
        radon_image; % the final FRT result
        radon_image_trans; % final transposed FRT
        
        % additional information/statistics:
        last_snr;
        snr_values; % a list of all SNR values found in this run (use "reset" to clear them)
        latest_thresh; % used in "purgeStreaks"
        
        subtracted_image; % image after PSF filter and after subtracting any found streaks (use_recursive to make sure all streaks are found)
        
    end
    
    properties(SetAccess=protected) % internal variables
        
        % lazy loaded from input_var when needed:
        noise_var; % scalar value: either given as scalar or the median value of the variance map
        radon_var_uni = {}; % if only a variance scalar is given, save a uniform variance map and multiply it by noise_var
        radon_var_map = {}; % if a 2D var-map is given, need to calculate an actual Radon var-map
        var_size; % size of image used for making the variance map/uniform
        var_was_expanded = 0; % keep track of the original size of var-map (if it was expanded)
        
        % lazy loaded from input_psf when needed: 
        psf; % psf map, either given as such or made using util.img.gaussian2
        psf_sigma; % scalar psf width (sigma parameter), either given or calculated using util.img.gaussfit. 
        
    end
    
    properties % switches/controls
        
        % image preprocessing
        use_subtract_mean = 1; % verify the mean of each image is zero
        use_subtract_median = 1; % verify the mean of each image is zero
        use_conv = 1; % use match filter on data (using the given PSF)        
        use_crop_image = 0; % crop image to fit powers of 2 (does not enlarge)
        crop_size = 2048; % what image size to crop to
        use_point_source_removal = 1;
        num_point_sources = 50; % how many point sources can we remove from the image
        
        % search options
        use_short = 1; % search for short streaks
        min_length = 32; % for finding short streaks (streak can be 1/cos(th) or 1/sin(th) larger)
        threshold = 10; % in units of S/N
        use_autothresh = 0; % find the threshold by comparing to SNR distribution
        use_recursive = 0; % search for multiple streaks (recursive search: remove line and do another FRT)
        recursion_depth = 10; % how many transforms are we willing to do to catch distinct lines (this is for each transposition)
        use_only_one = 1; % if you happen to get two streaks, one for transposed and one for untransposed, choose only the best one. Ignored when using recursion to find multiple streaks. 
        subtract_psf_widths = 3; % the width to subtract around the found position (when subtracting streaks). Units of PSF sigma
        autothresh_num_sigma = 5; % when doing sigma clipping and fitting the extreme-value distribution, how many "sigmas" above mean do you want to cut (used in "purgeStreaks")
        
        % post-Radon processing (i.e. exclusions)
        use_exclude = 1; % if you want to remove the row/column noise
        exclude_dx = [-1 1]*50; % how much to exclude from the transposed radon image
        exclude_dy = []; % how much to exclude from the non-transposed radon image
        
        % memory considerations
        use_save_images = 1; % use this to keep a copy of the original/final image in each streak
        use_clear_memory = 1; % use this to clear var maps and final images before saving object (can be recalculated using "recalculateFinalImages")
        
        % just to keep track of where the detection is made
        filename; % what is the filename from which we found the streak
        batch_num; % what batch the streak was found in (typically batch==file)
        frame_num; % which frame inside the batch
       
        % default display settings
        show_bit = 1; % if you want to display the results of each new input image
        display_index = 1; % what streak should be displayed
        display_which = 'current'; % can also be "previous" (streaks vs. prev_streaks)
        line_offset = 25; % default display line offset of the found streak in the input image
        rect_size = 250; % default display rectangle in the Radon image
        display_monochrome = 0; % display as monochrome, black-is-bright image
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
    end
    
    properties(Hidden=true)
        
        im_size_tr; % image size or transposed image size (depending on the last thing input to the FRT)
        
        psf_sigma_default = 2; % when no PSF is given, assume this is the width parameter "sigma"
        noise_var_default = 1; % when no variance is given, assume uniform variance map with value 1
        
        % these are filled at construction:
        default_min_length; 
        default_pixel_resolution;
        default_threshold;
        default_max_lines;
        default_crop_size;        
        default_exclude_dx;
        default_exclude_dy;
        
        version = 2.04;
        
    end
    
    methods % constructor
        
        function obj = Finder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'radon.Finder') % copy-constructor
                
                if obj.debug_bit, fprintf('Finder copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(varargin{1});
                
                if length(varargin)>1
                    obj.parse(varargin{2:end});
                end
                
            else

                if obj.debug_bit, fprintf('Finder constructor v%4.2f\n', obj.version); end

                util.oop.save_defaults(obj); % put default values for "X" into "default_X"
                
                obj.initialize;

                obj.parse(varargin{:}); % read any inputs from the constructor arguments

            end
            
        end
        
    end
    
    methods % reset/clear
        
        function initialize(obj) % should be done only when starting the survey/project/loadin matlab
            
            util.oop.load_defaults(obj); % load default value from "default_X" to "X"
            
            obj.timing_data = util.time.TimingData;
            
            obj.reset;
            
        end
        
        function reset(obj) % do this at the start of each data run

            obj.clearPSF;
            obj.clearVariance;
            
            obj.snr_values = [];
            obj.latest_thresh = [];

            obj.resetStreaks;
            
            obj.timing_data.reset;
            
        end
        
        function resetStreaks(obj) % clears streaks and prev_streaks
           
            obj.prev_streaks = radon.Streak.empty;
            obj.streaks = radon.Streak.empty;
                        
        end
        
        function clearPSF(obj) % get rid of the PSF image and scalar width
            
            obj.input_psf = [];
            obj.psf = [];
            obj.psf_sigma = [];
            
        end
        
        function clearVariance(obj) % get rid of all the variance maps and scalars
            
            obj.input_var = [];
            obj.var_size = [];
            obj.noise_var = [];
            obj.radon_var_uni = {};
            obj.radon_var_map = {};
            
            obj.timing_data.clear('make var');
            
        end
        
        function clear(obj) % clear the images and Radon images, and reset the timing data (call this for each batch)
            
            obj.input_images = [];
            obj.im_size = [];
            obj.im_size_tr = [];
            obj.radon_image = [];
            obj.radon_image_trans = [];
            obj.last_snr = [];
            obj.streaks = radon.Streak.empty;
            obj.last_streak = radon.Streak.empty;
            
            % we do not call timing_data.clear, because the variance map is sometimes only calculated per run (not per batch)
            obj.timing_data.clear('subtraction');
            obj.timing_data.clear('convolution');
            obj.timing_data.clear('show');
            obj.timing_data.clear('scan');
            obj.timing_data.clear('frt');
            
        end
        
    end
    
    methods % getters
        
        function val = get.im_size(obj) % verify the return value is a 2-element vector
            
            val = util.vec.imsize(obj.im_size);
            
        end
        
        function val = get.noise_var(obj) % always returns a scalar variance
            
            if isempty(obj.noise_var)
            
                if isempty(obj.input_var)
                    obj.noise_var = obj.noise_var_default;
                elseif isscalar(obj.input_var)                
                    obj.noise_var = obj.input_var;
                else
                    obj.noise_var = util.stat.median2(obj.input_var);
                end
                
            end
            
            val = obj.noise_var;
                        
        end
        
        function val = get.psf_sigma(obj) % always returns a scalar PSF width
            
            if isempty(obj.psf_sigma)

                if isempty(obj.input_psf)
                    obj.psf_sigma = obj.psf_sigma_default; 
                elseif isscalar(obj.input_psf)                
                    obj.psf_sigma = obj.input_psf;
                else
                    [p, res] = util.img.gaussfit(obj.input_psf, 'sigma_x', 0, 'amp', 0); % fit to a symmetric PSF with unknown normalization...
                    % should we test the residual first??
                    obj.psf_sigma = p(1); 
                end

            end
            
            val = obj.psf_sigma;
            
        end
        
        function val = get.psf(obj) % always returns a PSF image
            
            if isempty(obj.psf)
                
                if isempty(obj.input_psf) || isscalar(obj.input_psf)
                    obj.psf = util.img.gaussian2(obj.psf_sigma);
                    obj.psf = obj.psf./sqrt(util.stat.sum2(obj.psf.^2));
                else
                    obj.psf = obj.input_psf;
                    obj.psf = obj.psf./sqrt(util.stat.sum2(obj.psf.^2));
                end
                
            end
            
            val = obj.psf;

%             val = val./sqrt(sum(sum(val).^2)); % old version normalization 

        end
                
        function val = use_expand(obj) % do we need to expand the image before FRT (for short streaks: no)
            
            if obj.use_short
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = get.radon_var_uni(obj) % (lazy load) cell array of uniform variance maps of all partial transforms (for both transpositions)
            
            if isempty(obj.im_size)
                error('You are asking for a uniform var-map without giving an image size!');
            end
            
            if isempty(obj.radon_var_uni) && ( isempty(obj.input_var) || isscalar(obj.input_var) )
                
                obj.var_size = obj.im_size;
                
                obj.timing_data.start('make var');

                obj.radon_var_uni{1} = radon.frt(ones(obj.var_size), 'partial', 1, 'expand', obj.use_expand);
                if obj.im_size(1)==obj.im_size(2) % if the uniform map is square, the transposed element is the same as the un-transposed
                    obj.radon_var_uni{2} = obj.radon_var_uni{1}; 
                else
                    obj.radon_var_uni{2} = radon.frt(ones(obj.var_size), 'partial', 1, 'expand', obj.use_expand, 'transpose', 1);
                end
                
                obj.var_was_expanded = obj.use_expand; % keep track of whether the radon-var-map was expanded
                
                obj.timing_data.finish('make var');

            end
            
            val = obj.radon_var_uni;
            
        end
        
        function val = get.radon_var_map(obj) % (lazy load) cell array of non-uniform variance maps of all partial transforms (for both transpositions)
           
            if isempty(obj.radon_var_map) && ~isempty(obj.input_var) && ~isscalar(obj.input_var)
                
                obj.var_size = util.vec.imsize(obj.input_var);
                
                obj.timing_data.start('make var');

                obj.radon_var_map{1} = radon.frt(obj.input_var, 'partial', 1, 'expand', obj.use_expand);
                obj.radon_var_map{2} = radon.frt(obj.input_var, 'partial', 1, 'expand', obj.use_expand, 'transpose', 1);
                
                obj.var_was_expanded = obj.use_expand; % keep track of whether the radon-var-map was expanded
                
                obj.timing_data.finish('make var');

            end
            
            val = obj.radon_var_map;
            
        end
        
    end
    
    methods % setters
        
        function cycleDisplayWhich(obj)
                        
            import util.text.cs;
            
            if cs(obj.display_which, 'current')
                obj.display_which = 'prev';
            else
                obj.display_which = 'current';
            end
            
        end
        
        function set.input_images(obj, val) % also set the im_size parameter
            
            obj.input_images = val;
            
%             obj.im_size = util.vec.imsize(val);
            
        end
        
        function set.input_var(obj, val) % also clears the lazy loaded "noise_var" and "radon_var_map"
            
            if ~isempty(val) && ~isscalar(val)
                
                if ~isempty(obj.input_images) && (size(val,1)~=obj.im_size(1) || size(val,2)~=obj.im_size(2))
                    error('Size mismatch between image size %s and variance size %s', util.text.print_vec(obj.im_size, 'x'), util.text.print_vec(size(val), 'x'));
                end
                
                obj.radon_var_uni = {}; % if we are working with var-map we don't need this
                
            end
            
            obj.input_var = val;
            
            obj.noise_var = []; % will be lazy loaded
            obj.radon_var_map = {}; % will be lazy loaded
            
        end
        
        function set.input_psf(obj, val) % also clears the lazy loaded PSF variables
            
            if ~isempty(val)
                obj.clearPSF;
                obj.input_psf = val;
            end
            
        end
        
        function set.exclude_dx(obj, val) % give a scalar for symmetric bounds, or a 2-element vector
            
            if isscalar(val)
                obj.exclude_dx = [-abs(val) abs(val)];
            elseif length(val)==2
                obj.exclude_dx = val;
            elseif isempty(val)
                obj.exclude_dx = [];
            else
                error(['must input a scalar or two element vector to "exclude dx". ' util.text.printsize(val)]);
            end
            
        end
        
        function set.exclude_dy(obj, val) % give a scalar for symmetric bounds, or a 2-element vector
            
            if isscalar(val)
                obj.exclude_dy = [-abs(val) abs(val)];
            elseif length(val)==2
                obj.exclude_dy = val;
            elseif isempty(val)
                obj.exclude_dy = [];
            else
                error(['must input a scalar or two element vector to "exclude dy". ' util.text.printsize(val)]);
            end
            
        end
        
    end
        
    methods % getting the best streak/SNR/length
             
        function val = bestSNR(obj) % best SNR for all streaks in this frame...
            
            import util.stat.max2;
            
            if isempty(obj.streaks)
                if isempty(obj.radon_image)
                    val = 0;
                else
                    val = max(max2(obj.radon_image), max2(obj.radon_image_trans));
                end
            else
                val = max([obj.streaks.snr]);
            end

        end
        
        function val = bestLength(obj) % length of the best streak
            
            if isempty(obj.streaks)
                val = 0;
            else
                val = max([obj.streaks.L]);
            end
            
        end
        
        function val = bestStreak(obj) % streak with the highest SNR
            
            [~,idx] = max([obj.streaks.snr]);
            val = obj.streaks(idx);
            
        end
        
        function val = bestPrevStreak(obj) % streak with the highest SNR from all streaks seen in this run
            
            [~,idx] = max([obj.prev_streaks.snr]);
            val = obj.streaks(idx);
            
        end
        
    end
        
    methods % calculations
    
        function parse(obj, varargin) % reconize some keyword-value pairs as inputs and parameters
                      
            import util.text.cs;
            import util.text.parse_bool;
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'variance')
                    obj.input_var = val;
                elseif cs(key, 'psf')
                    obj.input_psf = val;
                elseif cs(key, {'short', 'use_short'})
                    obj.use_short = parse_bool(val);                    
                elseif cs(key, {'min_length', 'minimal_length', 'length'})
                    obj.min_length = val;
                elseif cs(key, 'threshold')
                    obj.threshold = val;
                elseif cs(key, {'use_recursive', 'recursive'})
                    obj.use_recursive = parse_bool(val);
                elseif cs(key, 'depth_recursion')
                    obj.recursion_depth = val;
                elseif cs(key, {'exclude', 'use_exclude'})
                    obj.use_exclude = parse_bool(val);
                elseif cs(key, {'dx'})
                    obj.exclude_dx = val;
                elseif cs(key, {'dy'})
                    obj.exclude_dy = val;
                end
                
            end
            
        end
        
        function V = getRadonVariance(obj, transpose, step) % gets the correct radon var map for this log-step and transposition (creates a var-map if needed)
           
            import util.stat.sum2;
            
            if nargin<2 || isempty(transpose)
                transpose = 0;
            end
            
            if nargin<3 || isempty(step)
                if transpose
                    step = ceil(log2(obj.im_size(1)));
                else
                    step = ceil(log2(obj.im_size(2)));
                end
            end
            
            if ~isempty(obj.var_size) && (obj.im_size(1)~=obj.var_size(1) || obj.im_size(2)~=obj.var_size(2) || obj.use_expand~=obj.var_was_expanded)
                % if we changed the size/expansion of the image
                obj.timing_data.clear('make var');
                obj.radon_var_uni = {};
                obj.radon_var_map = {};
            end
            
            if isempty(obj.input_var) || isscalar(obj.input_var)
                V = obj.radon_var_uni{transpose+1}{step}*obj.noise_var; % lazy load through getter function
            else
                V = obj.radon_var_map{transpose+1}{step}; % lazy load through getter function
            end
            
        end
        
        function input(obj, varargin) % THIS IS THE USEFUL METHOD (input images, variance and PSF to find streaks)
            
            import util.text.cs;
            import util.img.crop2size;
            
            input = util.text.InputVars; % parses varargin-pairs
            input.input_var('images', []);
            input.input_var('variance', [], 'noise_variance');
            input.input_var('psf', [], 'point spread function');
            input.input_var('filename', []);
            input.input_var('batch_number', []);
            input.use_ordered_numeric = 1;
            input.scan_vars(varargin{:});
            
            if isempty(input.images) % don't clear/continue the function without images...
                return;
            end
            
            obj.clear;
            
            if obj.use_crop_image % use this to crop down to power of 2
                obj.input_images = crop2size(input.images, obj.crop_size);
            else
                obj.input_images = input.images;
            end
            
            obj.im_size = util.vec.imsize(obj.input_images);
            
            if ~isempty(input.variance) % if no variance is given, continue with what you have (or default value)
                if obj.use_crop_image
                    obj.input_var = crop2size(input.variance, obj.crop_size);    
                else
                    obj.input_var = input.variance;
                end
            end
                
            if ~isempty(input.psf) % if no PSF is given, continue with what you have (or default value)
                obj.input_psf = input.psf;
            end
            
            % housekeeping variables
            obj.filename = input.filename;
            obj.batch_num = input.batch_number;
            
            % need these to normalize the radon image
            V = obj.getRadonVariance(0); % if radon_var_map exists, use it. otherwise, get a radon_var_uni instead
            VT = obj.getRadonVariance(1);
            
            V = permute(V, [1,3,2]);
            VT = permute(VT, [1,3,2]);
            
            N = pow2(nextpow2(obj.im_size(1)));
            th = atand((-N+1:N-1)./obj.im_size(1)); % angle range, in degrees
            G = max(abs(cosd(th)), abs(sind(th))); % geometric factor correction to the S/N
            NT = pow2(nextpow2(obj.im_size(2)));
            thT = atand((-NT+1:NT-1)./obj.im_size(2)); % angle range, in degrees
            GT = max(abs(cosd(thT)), abs(sind(thT))); % geometric factor correction to the S/N for the transposed image
            
            for ii = 1:size(obj.input_images,3) % loop on all images in the batch
                
                obj.frame_num = ii;
                obj.last_snr = [];
                
                if obj.debug_bit
                    fprintf('running streak detection on batch= % 2d | frame= % 3d... ', obj.batch_num, obj.frame_num);
                end

                I = obj.preprocess(obj.input_images(:,:,ii), V, VT, G, GT); % subtract background, remove point sources, get rid of NaN values
                
                if isscalar(I)
                    disp(['adapting the SNR value ' num2str(I) ' from 1st pass']);
                    new_snr = I; % this only happens if a point source search triggers a high threshold FRT search and finds a really bright streak... 
                else
                    new_snr = obj.sendToFRT(obj.input_images(:,:,ii), I, V, VT, G, GT); % send to FRT, transpose, remove streaks, find best SNR
                end
                
                if obj.debug_bit
                    fprintf(' | length= %d | bestSNR= %f\n', round(obj.bestLength), new_snr);
                end
                
                obj.snr_values = [obj.snr_values new_snr];
                
                if ~isempty(obj.streaks)
                    if obj.debug_bit, disp(['saving streaks with SNR= ' util.text.print_vec(new_snr)]); end
                    obj.prev_streaks = [obj.prev_streaks obj.streaks];                    
                end
                
                obj.timing_data.start('show');
                
                if obj.show_bit && obj.gui.check
                    obj.show;
                end
                
                obj.timing_data.finish('show');
                
                
                drawnow;
                
            end
            
        end
        
        function I_filtered = preprocess(obj, I_original, V, VT, G, GT)

            import util.vec.pick_index;
            import util.stat.max2;
            import util.stat.mean2;
            import util.stat.median2;
            
            obj.timing_data.start('subtraction');

            I = I_original;
                        
            if obj.use_subtract_median
                I = I - median2(I);
            end
            
            if obj.use_subtract_mean
                I = I - mean2(I);
            end

            I(isnan(I)) = 0; % If you cut the stars out, replace them with NaNs, then subtract mean (see above) then replace NaNs with zeros.

            obj.timing_data.finish('subtraction');

            obj.timing_data.start('convolution');

            if obj.use_conv && ~isempty(obj.psf)
                I_filtered = filter2(pick_index(obj.psf, obj.frame_num, 3), I);
            else
                I_filtered = I;
            end

            obj.timing_data.finish('convolution');

            obj.timing_data.start('point_search');

            if obj.use_point_source_removal

                I_filtered_var = filter2(pick_index(obj.psf, obj.frame_num, 3), I./sqrt(obj.input_var));
                
                for jj = 1:obj.num_point_sources

                    if mod(jj,10)==0, fprintf('% 4d\n', jj); end
                    
                    [mx, idx] = max2(I_filtered_var);

                    if mx>obj.threshold
                        
                        if jj==1 % maybe the "point sources" are actually pieces of a very bright streak, so try to run high-threshold detection on them
                            if obj.debug_bit, fprintf(' (PS peak: %f... triggered 1st pass with high thresh) ', mx); end
                            new_snr = obj.sendToFRT(I_original,I_filtered,V,VT,G,GT);
                            if new_snr>mx.*sqrt(obj.min_length./obj.psf_sigma) % the streak has to be brighter than the point source! 
                                I = new_snr;
                                return; % this happens only if we found a very bright streak that triggered the point source detector
                            end
                        end
                        
                        delta = ceil(obj.psf_sigma.*2);
                        x_low = idx(2)-delta;
                        if x_low<1, x_low = 1; end
                        x_high = idx(2)+delta;
                        if x_high>size(I,2), x_high=size(I,2); end
                        
                        y_low = idx(1)-delta;
                        if y_low<1, y_low = 1; end
                        y_high = idx(1)+delta;
                        if y_high>size(I,1), y_high=size(I,1); end
                        
                        I(y_low:y_high,x_low:x_high) = NaN;
                        
                        util.plot.show(I, 'dyn', 100); drawnow; % debug only! 
                        
                    else
                        break; % if there are no more bright point sources, just move on
                    end

                end

                if jj>1 % if any sources were removed refind the mean and then get rid of the NaNs

                    if obj.use_subtract_median
                        I = I - median2(I);
                    end
                    
                    if obj.use_subtract_mean
                        I = I - mean2(I);
                    end

                    I(isnan(I)) = 0; % If you cut the stars out, replace them with NaNs, then subtract mean (see above) then replace NaNs with zeros.

                    if obj.use_conv && ~isempty(obj.psf) % also update the filtered image
                        I_filtered = filter2(pick_index(obj.psf, obj.frame_num, 3), I);
                    else
                        I_filtered = I;
                    end

                end
                
            end

            obj.timing_data.finish('point_search');

        end
        
        function new_snr = sendToFRT(obj, I_original, I_filtered, V, VT, G, GT)

            import util.stat.sum2;
            
            obj.timing_data.start('frt');

            obj.streaks = radon.Streak.empty;
            obj.last_streak = radon.Streak.empty;
            R = radon.frt(I_filtered, 'finder', obj, 'expand', obj.use_expand);
            temp_streaks = obj.streaks; % keep this here while filling the list with the transposed image
            obj.streaks = radon.Streak.empty;
            obj.last_streak = radon.Streak.empty;

            if obj.use_only_one==0 || obj.use_recursive
                RT = radon.frt(obj.subtracted_image, 'finder', obj, 'expand', obj.use_expand, 'trans', 1);
            else
                RT = radon.frt(I_filtered, 'finder', obj, 'expand', obj.use_expand, 'trans', 1);
            end

            obj.streaks = [temp_streaks obj.streaks]; % combine the two lists

            if obj.use_only_one && ~obj.use_recursive % if you get two streaks (one for each transpose) but only want the best one
                obj.streaks = obj.bestStreak;
            end

            obj.timing_data.finish('frt');

            obj.radon_image = R./sqrt(V.*sqrt(sum2(obj.psf.^2).*G).*sum2(obj.psf)); % this normalization makes sure the result is in units of S/N, regardless of the normalization of the PSF. 
            obj.radon_image_trans = RT./sqrt(VT.*sqrt(sum2(obj.psf.^2).*GT).*sum2(obj.psf)).*GT; % one makes sure the PSF is unity normalized, the other divides by the amount of noise "under" the PSF

            % check if we want to save a copy of the subframe and input/final image in each of the streaks... 
            if obj.use_save_images
                for jj = 1:length(obj.streaks)

                    obj.streaks(jj).input_image = I_original;

                    if obj.streaks(jj).transposed
                        obj.streaks(jj).radon_image = obj.radon_image_trans;
                    else
                        obj.streaks(jj).radon_image = obj.radon_image;
                    end

                end
            else
                for jj = 1:length(obj.streaks)
                    obj.streaks(jj).input_image = [];
                    obj.streaks(jj).subframe = [];
                    obj.streaks(jj).radon_image = [];
                end
            end

            new_snr = max(obj.last_snr, obj.bestSNR);

        end
        
        function scan(obj, M, transpose) % do this on the final Radon image and on partial Radon matrices (for short streaks)
        % scans the given subsection of the FRT for streaks. 
        % M is the FRT matrix (can be full-FRT of subsection)
        % m is how many binary addition levels were performed to get it,
        % where we calculate m from the number of columns we have. 
            
            import util.stat.sum2;            
            import util.stat.maxnd;
            
            obj.timing_data.start('scan');
                        
            if nargin<3 || isempty(transpose)
                transpose = 0;
            end
                        
            if size(M,3)==1 % if we are given the matrix after it was permuted from 3D to 2D
                M = permute(M, [1,3,2]);
            end
            
            m = log2(size(M,3)+1)-1; % what logarithmic step we are in the FRT
            
            if obj.use_short==0 && 2^m<obj.im_size_tr(2) % if you don't want to scan partial matrices / short streaks
                return; 
            end
            
            if 2^m<obj.min_length/8 % don't bother with short lines (take three step below the given min_length, to eliminate point sources)
                return;
            end
            
            V = obj.getRadonVariance(transpose, m);
%             th = atand((-obj.im_size_tr(1)+1:obj.im_size_tr(1)-1)./obj.im_size_tr(1)); % angle range, in degrees
            S = (size(M,3)+1)/2;
            th = atand((-S+1:S-1)./S);
            G = max(abs(cosd(th)), abs(sind(th))); % geometric factor correction to the S/N
            G = util.vec.topages(G);
            
            SNR = double(M)./sqrt(V*sqrt(sum2(obj.psf.^2))*sum2(obj.psf).*G); % this should be normalized to SNR units... (see the same normalization in "input")
            SNR_final = SNR; % copy not to be saved (e.g. exclusions in the middle)
            
            if obj.use_exclude % get rid of vertical/horizontal lines
                
                if transpose && ~isempty(obj.exclude_dx)
                    offset = (size(M,3)+1)/2; % index of dy=0, also how many pixels are for angles 0<=th<=45 in this subframe
                    scale = obj.im_size_tr(2)/offset; % scaling factor for smaller subframes
                    idx1 = offset + ceil(obj.exclude_dx(1)/scale);
                    idx2 = offset + floor(obj.exclude_dx(2)/scale);
                    SNR_final(:, :, idx1:idx2) = 0;
                elseif ~isempty(obj.exclude_dy)
                    offset = (size(M,3)+1)/2;
                    scale = obj.im_size_tr(1)/offset;
                    idx1 = offset + ceil(obj.exclude_dy(1)/scale);
                    idx2 = offset + floor(obj.exclude_dy(2)/scale);
                    SNR_final(:, :, idx1:idx2) = 0;
                end
                                
            end
            
            [mx, idx] = maxnd(SNR_final);
            
            if isempty(obj.last_snr) || mx>obj.last_snr
                obj.last_snr = mx;
            end
            
            if mx>obj.threshold && (isempty(obj.last_streak) || mx>obj.last_streak.snr) % if a shorter streak is already saved, only update it if the new streak has higher SNR
            
                obj.last_streak = radon.Streak(obj, SNR, m, transpose, M(idx(1), idx(2), idx(3)), idx);
                
            end
            
            obj.timing_data.finish('scan');
                      
        end
        
        function finalizeFRT(obj, M_in, transpose, radon_image) % call this at the end of the FRT function, if there is a finder
            
            obj.subtracted_image = M_in; % must have something in this image
            
            if ~isempty(obj.last_streak) && obj.last_streak.radon_dx < obj.min_length % kill streaks that are too short (e.g. bright point-sources)
                obj.last_streak = radon.Streak.empty;
%                 obj.last_snr = 0;
            end

            if ~isempty(obj.last_streak) % if we found a streak
                
                obj.last_streak.input_image = M_in;
                obj.radon_image = radon_image;
                obj.streaks(end+1) = obj.last_streak;
                
                if ~obj.use_only_one
                    obj.subtracted_image = obj.streaks(end).subtractStreak(M_in, obj.subtract_psf_widths); % if any streaks are found, subtract them first... 
                end
                
                if obj.use_recursive && length(obj.streaks)<obj.recursion_depth % keep calling FRT until all streaks are found...
                    
                    if obj.debug_bit>9 % use debug_bit=10 to watch the streaks getting removed (make sure no edges are left!)
                        util.plot.show(obj.subtracted_image, 'bias', 0, 'dyn', 100);
                        title(sprintf('num streaks found: %d | transpose: %d', length(obj.streaks), transpose));
                        drawnow;
                        pause(1);
                    end
                    
                    radon.frt(obj.subtracted_image, 'transpose', transpose, 'expand', obj.use_expand, 'finder', obj); % recursively call FRT again with the same input image, without the found streak... 
                    
                end
                
            end

        end
                
        function purgeStreaks(obj) % call this at the end of each run to find an optimal threshold and then get rid of all streaks below it (from "prev_streaks" too)
            
            snr = obj.snr_values(obj.snr_values>0 & obj.snr_values<100); % kill zero values (too short) and very bright values
            
            [mu, sigma] = util.stat.sigma_clipping(snr, 'dist', 'max'); % iteratively remove outliers from an extreme-value distribution
            
            obj.latest_thresh = mu+obj.autothresh_num_sigma*sigma; % usually autothresh_num_sigma=5
            
            obj.prev_streaks = obj.prev_streaks([obj.prev_streaks.snr]>obj.latest_thresh); % cut out the streaks with the new threshold
            
            for ii = 1:length(obj.prev_streaks)
                if obj.prev_streaks(ii).threshold<obj.latest_thresh
                    obj.prev_streaks(ii).threshold = obj.latest_thresh;
                end
            end
            
        end

    end
        
    methods % plotting tools / GUI
        
        function makeGUI(obj) % fills the GUI object if needed
           
            if isempty(obj.gui)
                obj.gui = radon.gui.FinderGUI(obj);
            end
            
            obj.gui.makeGUI;
            
        end
        
        function show(obj, varargin) % displays the input/Radon image and some streak parameters (usually for the GUI)
            
            import util.text.cs;
            
            % parse varargin pairs:
            
            if cs(obj.display_which, 'current')
                if isempty(obj.streaks)
                    return;
                end
                s = obj.streaks(obj.getStreakIndex);
            elseif cs(obj.display_which, 'previous')
                if isempty(obj.prev_streaks)
                    return;
                end
                s = obj.prev_streaks(obj.getStreakIndex);
            else
                error(['Unknown display_which option "' obj.display_which '", use current/previous...']);
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                parent = obj.gui.panel_image;
            else
                parent = [];
            end
            
            s.show('parent', parent, 'line_offset', obj.line_offset, 'rect_size', obj.rect_size, 'monochrome', obj.display_monochrome);
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
        end
        
        function showAllStreaks(obj, varargin)

            import util.text.cs;

            if cs(obj.display_which, 'current')
                N = length(obj.streaks);
            elseif cs(obj.display_which, 'previous')
                N = length(obj.prev_streaks);
            else
                error(['Unknown display_which option "' obj.display_which '", use current/previous...']);
            end
            
            if N==0
                return;
            end
            
            obj.display_index = 1;
            
            for ii = 1:N
                obj.show(varargin{:});
                pause(1);
                obj.nextDisplayIndex;
            end
            
        end
        
        function val = getStreakIndex(obj) % finds the index of the best streak (for "streaks" or "prev_streaks")
            
            import util.text.*;
            
            if cs(obj.display_which, 'current')
                
                if isempty(obj.display_index)
                    [~,val] = max([obj.streaks.snr]);
                elseif obj.display_index>length(obj.streaks) || obj.display_index<1
                    val = 1;
                else
                    val = obj.display_index;
                end
                
            elseif cs(obj.display_which, 'previous')
               
                if isempty(obj.display_index)
                    [~,val] = max([obj.prev_streaks.snr]);
                elseif obj.display_index>length(obj.prev_streaks) || obj.display_index<1
                    val = 1;
                else
                    val = obj.display_index;
                end
                
            else
                error(['Unknown display_which option "' obj.display_which '", use current/previous...']);
            end
            
        end
        
        function resetDisplayIndex(obj) % will be lazy loaded by "getStreakIndex"
            obj.display_index = [];
        end
        
        function prevDisplayIndex(obj) % push index back
            
            import util.text.*;
            
            obj.display_index = obj.getStreakIndex - 1;
            
            if obj.display_index<1

                if cs(obj.display_which, 'current')                
                    obj.display_index = length(obj.streaks);
                elseif cs(obj.display_which, 'previous')
                    obj.display_index = length(obj.prev_streaks);
                else
                    error(['Unknown display_which option "' obj.display_which '", use current/previous...']);
                end

            end
            
        end
        
        function nextDisplayIndex(obj) % push index forward
            
            import util.text.*;
            
            obj.display_index = obj.getStreakIndex + 1;
            
            if cs(obj.display_which, 'current')                
                N = length(obj.streaks);
            elseif cs(obj.display_which, 'previous')
                N = length(obj.prev_streaks);
            else
                error(['Unknown display_which option "' obj.display_which '", use current/previous...']);
            end
            
            if obj.display_index>N
                obj.display_index = 1;
            end
            
        end
    end    
    
    methods % utilities
       
        function saveas(obj, filename, varname) % saves the finder as MAT file
                       
            if nargin<3 || isempty(varname)
                [~, varname, ext] = fileparts(filename);
            else
                [~, ~, ext] = fileparts(filename);
            end
            
            if isempty(ext)
                filename = [filename '.mat'];
            end
            
            if obj.debug_bit
                disp(['saving Finder object "' varname '" in file: ' filename]);
            end
            
            if obj.use_clear_memory
                obj.clearMemory;
            end
            
            save_struct.(varname) = obj;
            
            save(filename, '-struct', 'save_struct', '-v7.3');
            
        end
        
        function clearMemory(obj) % gets rid of images/Radon images/var maps to save memory (e.g. before saving to disk)

            if obj.debug_bit, disp('reducing memory by clearing final images...'); end

            obj.input_images = [];
            obj.radon_image = [];
            obj.radon_image_trans = [];
            obj.radon_var_uni = {};
            obj.radon_var_map = {};

            for ii = 1:length(obj.streaks)
                obj.streaks(ii).radon_image = [];
                obj.streaks(ii).subframe = [];
            end

            for ii = 1:length(obj.prev_streaks)
                obj.prev_streaks(ii).radon_image = [];
                obj.prev_streaks(ii).subframe = [];
            end
            
        end
        
        function recalculateFinalImages(obj) % rebuild all that was erased by "clearMemory"
            
            import util.stat.sum2;
            
            V = permute(obj.getRadonVariance(0), [1,3,2]);
            VT = permute(obj.getRadonVariance(1), [1,3,2]);
            
            s_vec = [obj.streaks obj.prev_streaks]; % go over all streaks
            
            for ii = 1:length(s_vec)
                
                s = s_vec(ii);
                
                I = s.input_image;
                
                I(isnan(I)) = 0;
                
                if s.was_convolved
                    I = filter2(s.psf, I);
                end
                
                R_partial = radon.frt(I, 'expand', s.was_expanded, 'transpose', s.transposed, 'partial', 1);
                
                R = radon.partial2final(R_partial{end});
                
                if s.transposed
                    s.radon_image = R./sqrt(VT.*sum2(s.psf).*sqrt(sum2(s.psf.^2)));
                else
                    s.radon_image = R./sqrt(V.*sum2(s.psf).*sqrt(sum2(s.psf.^2)));
                end
                
                s.subframe = R_partial{s.radon_step}./sqrt(obj.getRadonVariance(s.transposed, s.radon_step));
                
            end
            
        end
        
    end
    
    methods(Static=true)
       
        function [streaks, snr_values] = input_static(finder, I, worker_number) % for using a parpool to divide many images to several workers
           
            tic;
            finder.input(I); 
            streaks = finder.prev_streaks;
            snr_values = finder.snr_values;
            disp(['worker: ' num2str(worker_number) ' is finished after ' num2str(toc)]); 
            
        end
        
    end
    
end

