classdef Finder < handle
% Streak finding tool. 
% This object can be used in two ways:
%   (1) Give it an image, using obj.input(image, ...)
%       This will do all the internal calculations needed to find streaks.
%       Optional arguments to "input" are:
%           -variance: give a scalar (average) variance *or* a variance map.
%                      If no variance is given, Finder will use the previous
%                      value/map or the default which is variance=1.
%           -psf: give the point spread function as a scalar width (in this 
%                 case the PSF is a 2D Gaussian with sigma equal to the input 
%                 scalar) *or* a map for the image. 
%                 If no PSF is given, Finder will use the previous
%                 value/map or the default which is psf=1.
%           -filename: housekeeping parameter. Helps keep track of where 
%                      each streak was found.
%           -batch number: housekeeping parameter. Helps keep track of where 
%                          each streak was found.
%           -frame number: housekeeping parameter. Helps keep track of where 
%                          each streak was found.
%           -section number: housekeeping parameter. Helps keep track of where 
%                            each streak was found (in case you split images 
%                            to sections).
%           -offset: the top-left corner of the section inside the larger 
%                    image. Input y-then-x of first point. Default is [1 1].
%
%   NOTES: *Finder will run FRT twice, once on each transpose of the input. 
%           By default it will only save the highest SNR streak it finds. 
%           To save both, set "use_only_one=0". 
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
%    NOTES: *In this way the FRT is done only on one transposition (choose 
%            which one using the frt(..., transpose=true/false) parameter). 
% 
% SWITCHES AND SETTINGS (see comments in properties block for more details)
%   -Image pre-processing: use_subtract_mean, subtract_median, use_conv, 
%                          use_crop, crop_size... 
%   -Search options: use_short, min_length, threshold. 
%   -Post processing: use_exclude, exclude_dx, exclude_dy
%   -memory considerations: use_save_images, use_clear_memory.
%   -housekeeping: filename, batch_num, frame_num, section_num.
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
        
        streak@radon.Streak; % latest streak found by a call to input(...) or by sending Finder as input to FRT. If "use_only_one==0" then can contain 2 streaks. 
        prev_streaks@radon.Streak; % a vector of all streaks since the last call to "reset" (all streaks per run)
                        
    end
        
    properties % inputs/outputs
                
        % images:
        image; % what is given as input to the FRT.         
        im_size; % can be input as scalar, will be output as two-element vector.
        
        original_image; % given by user for displaying the streaks
        offset_original; % y-then-x offset of top-left corner of the section/crop relative to original image. 
        was_convolved = 0; % if the input image is convolved before it was given to Finder.
        
        % input var/psf:
        input_var; % original variance map/scalar 
        input_psf; % original PSF as given (or width parameter "sigma" of PSF)        
        
        % outputs:
        radon_image; % the final FRT result
        radon_image_trans; % final transposed FRT
        
        radon_var; % the variance map, in radon space
        radon_var_trans; % the transposed variance map, in radon space
        
        % additional information/statistics:
        last_snr; % the latest streak's SNR
        snr_values; % a list of all SNR values found in this run (use "reset" to clear them). Includes streaks that were not saved! 
                
        subtracted_image; % image after PSF filter and after subtracting any found streaks (use_recursive to make sure all streaks are found) <---- Do we need this anymore??
        
    end
    
    properties(SetAccess=protected) % internal variables
        
        % lazy loaded from input_var when needed:
        noise_var; % scalar value: either given as scalar or the median value of the variance map
        radon_var_uni = {}; % if only a variance scalar is given, save a uniform variance map and multiply it by noise_var
        radon_var_map = {}; % if a 2D var-map is given, need to calculate an actual Radon var-map
        var_size; % size of image used for making the variance map/uniform
        var_was_expanded = 0; % keep track of the original size of var-map (if it was expanded)
        
        % lazy loaded from input_var when needed:
        var_map; % map of the variance, same size as image_processed
        var_scalar; % either a scalar given to input_var, or the median of the var_map
        
        % lazy loaded from input_psf when needed: 
        psf; % psf map, either given as such or made using util.img.gaussian2
        psf_sigma; % scalar psf width (sigma parameter), either given or calculated using util.img.gaussfit. 
        
    end
    
    properties % switches/controls
        
        % image preprocessing
        use_preprocess = 1; % apply additional processing before image goes to FRT
        use_subtract_mean = 1; % verify the mean of each image is zero
        use_subtract_median = 0; % verify the median of each image is zero (supercedes use_subtract_mean)
        use_conv = 1; % use match filter on data (using the given PSF)
        use_crop = 0; % crop image to fit powers of 2 (does not enlarge)
        crop_size = 2048; % what image size to crop to
        
        % search options
        use_short = 1; % search for short streaks
        min_length = 32; % for finding short streaks (streak can be 1/cos(th) or 1/sin(th) larger)
        threshold = 10; % in units of S/N
        use_only_one = 1; % if you happen to get two streaks, one for transposed and one for untransposed, choose only the best one. Ignored when using recursion to find multiple streaks. 
        subtract_psf_widths = 5; % the width to subtract around the found position (when subtracting streaks). Units of PSF sigma <---- Do we need this anymore??
        
        % post-Radon processing (i.e. exclusions)
        use_exclude = 0; % if you want to remove the row/column noise
        exclude_dx = [-1 1]*50; % how much to exclude from the transposed radon image
        exclude_dy = []; % how much to exclude from the non-transposed radon image
        
        % memory considerations
        use_save_images = 1; % use this to keep a copy of the original/final image in each streak
        use_clear_memory = 1; % use this to clear var maps and final images before saving object (can be recalculated using "recalculateFinalImages")
        
        % just to keep track of where the detection is made
        filename; % what is the filename from which we found the streak
        batch_num; % what batch the streak was found in (typically batch==file)
        frame_num; % which frame inside the batch
        section_num; % which section inside the image was used (for sectioned images)
        
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
        
        psf_sigma_default = 1; % when no PSF is given, assume this is the width parameter "sigma"
        var_scalar_default = 1; % when no variance is given, assume uniform variance map with value 1
                
        % these are filled at construction:
        default_min_length; 
        default_pixel_resolution;
        default_threshold;
        default_crop_size;        
        default_exclude_dx;
        default_exclude_dy;
        
        version = 3.00;
        
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
            
            obj.reset;
            
        end
        
        function reset(obj) % do this at the start of each data run

            obj.clearPSF;
            obj.clearVariance;
            
            obj.snr_values = [];
            obj.resetStreaks;
            
        end
        
        function resetStreaks(obj) % clears streaks and prev_streaks
           
            obj.prev_streaks = radon.Streak.empty;
            obj.streak = radon.Streak.empty;
                        
        end
        
        function clearPSF(obj) % get rid of the PSF image and scalar width
            
            obj.input_psf = [];
            obj.psf = [];
            obj.psf_sigma = [];
            
        end
        
        function clearVariance(obj) % get rid of all the variance maps and scalars
            
            obj.input_var = [];
            obj.var_size = [];
            obj.var_map = [];
            obj.var_scalar = [];
            obj.noise_var = []; % redundant with var_scalar
            obj.radon_var_uni = {};
            obj.radon_var_map = {};
            obj.radon_var = [];
            obj.radon_var_trans = [];
                        
        end
        
        function clear(obj) % clear the images and Radon images, and reset the timing data (call this for each batch)
            
            obj.image = [];            
            obj.im_size = [];
            obj.im_size_tr = [];
            
            obj.original_image = [];
            obj.offset_original = [0 0];
            obj.was_convolved = 0;
            
            obj.radon_image = [];
            obj.radon_image_trans = [];
            
            obj.last_snr = 0;
            obj.streak = radon.Streak.empty;
            
        end
        
    end
    
    methods % getters
        
        function val = get.im_size(obj) % verify the return value is a 2-element vector
            
            val = util.vec.imsize(obj.im_size);
            
        end
        
        function val = get.noise_var(obj) % always returns a scalar variance
            
            if isempty(obj.noise_var)
            
                if isempty(obj.input_var)
                    obj.noise_var = obj.var_scalar_default;
                elseif isscalar(obj.input_var)                
                    obj.noise_var = obj.input_var;
                else
                    obj.noise_var = util.stat.median2(obj.input_var);
                end
                
            end
            
            val = obj.noise_var;
                        
        end
        
        function val = get.var_scalar(obj) % always returns a scalar variance
            
            if isempty(obj.var_scalar)
            
                if isempty(obj.input_var)
                    obj.var_scalar = obj.var_scalar_default;
                elseif isscalar(obj.input_var)
                    obj.var_scalar = obj.input_var;
                else
                    obj.var_scalar = util.stat.median2(obj.var_map);
                end
                
            end
            
            val = obj.noise_var;
                        
        end
        
        function val = get.var_map(obj)
            
            if isempty(obj.var_map)
                if isempty(obj.input_var) || isscalar(obj.input_var)
                    if ~isempty(obj.image)
                        obj.var_map = ones(size(obj.image)).*obj.var_scalar;
                    elseif ~isempty(obj.image)
                        obj.var_map = ones(size(obj.image)).*obj.var_scalar;
                    else
                        obj.var_map = [];
                    end
                else
                    obj.var_map = obj.input_var;
                end
            end
            
            val = obj.var_map;
            
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
                
                obj.radon_var_uni{1} = radon.frt(ones(obj.var_size), 'partial', 1, 'expand', obj.use_expand);
                if obj.im_size(1)==obj.im_size(2) % if the uniform map is square, the transposed element is the same as the un-transposed
                    obj.radon_var_uni{2} = obj.radon_var_uni{1}; 
                else
                    obj.radon_var_uni{2} = radon.frt(ones(obj.var_size), 'partial', 1, 'expand', obj.use_expand, 'transpose', 1);
                end
                
                obj.var_was_expanded = obj.use_expand; % keep track of whether the radon-var-map was expanded
                
            end
            
            val = obj.radon_var_uni;
            
        end
        
        function val = get.radon_var_map(obj) % (lazy load) cell array of non-uniform variance maps of all partial transforms (for both transpositions)
           
            if isempty(obj.radon_var_map) && ~isempty(obj.input_var) && ~isscalar(obj.input_var)
                
                obj.var_size = util.vec.imsize(obj.input_var);
                
                obj.radon_var_map{1} = radon.frt(obj.input_var, 'partial', 1, 'expand', obj.use_expand);
                obj.radon_var_map{2} = radon.frt(obj.input_var, 'partial', 1, 'expand', obj.use_expand, 'transpose', 1);
                
                obj.var_was_expanded = obj.use_expand; % keep track of whether the radon-var-map was expanded
                
            end
            
            val = obj.radon_var_map;
            
        end
        
        function val = getNormFactorPSF(obj)
            
            import util.stat.sum2;
            
            if (obj.use_conv && obj.use_preprocess) || obj.was_convolved
                val = sqrt(sum2(obj.psf.^2)).*sum2(obj.psf);
            else
                val = 1;
            end
        end
        
        function val = getGeometricFactor(obj, m_step)
        
%             S = 2^(m_step+1)-1;
            S = 2^m_step;
            th = atand((-S+1:S-1)./S);
            val = max(abs(cosd(th)), abs(sind(th))); % geometric factor correction to the S/N
            
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
        
        function set.input_var(obj, val) % also clears the lazy loaded "noise_var" and "radon_var_map"
            
            if ~isempty(val) && ~isscalar(val)
                
                if ~isempty(obj.image) && (size(val,1)~=obj.im_size(1) || size(val,2)~=obj.im_size(2))
                    error('Size mismatch between image size %s and variance size %s', util.text.print_vec(obj.im_size, 'x'), util.text.print_vec(size(val), 'x'));
                end
                
                obj.radon_var_uni = {}; % if we are working with var-map we don't need this
                
            end
            
            if isempty(val) || isempty(obj.input_var) || val~=obj.input_var
                
                obj.noise_var = []; % will be lazy loaded
                obj.var_scalar = []; % will be lazy loaded
                obj.var_map = []; % will be lazy loaded
                obj.radon_var_map = {}; % will be lazy loaded

                obj.input_var = val;
                
            end
            
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
             
        function val = bestSNR(obj) % best SNR, either from one of the streaks or the peak of the final radon image (or the transposed radon image)
            
            import util.stat.max2;
            
            if isempty(obj.streak) % no streaks: find the peak of the radon image / trans. radon image. 
                
                if ~isempty(obj.radon_image) 
                    R_max = max2(obj.radon_image);
                else
                    R_max = 0;
                end
                
                if ~isempty(obj.radon_image_trans)
                    RT_max = max2(obj.radon_image_trans);
                else
                    RT_max = 0;
                end
                
                val = max(R_max, RT_max);
                
            else
                val = max([obj.streak.snr]); % streaks may have a higher SNR if they are short! 
            end

        end
        
        function val = bestLength(obj) % length of the best streak (if none are found, set to zero...) 
            
            if isempty(obj.streak)
                val = 0;
            else
                val = max([obj.streak.L]);
            end
            
        end
        
        function val = bestStreak(obj) % streak with the highest SNR
            
            [~,idx] = max([obj.streak.snr]);
            val = obj.streak(idx);
            
        end
        
        function val = bestPrevStreak(obj) % streak with the highest SNR from all streaks seen in this run
            
            [~,idx] = max([obj.prev_streaks.snr]);
            val = obj.prev_streaks(idx);
            
        end
        
    end
        
    methods % calculations
    
        function streaks = input(obj, varargin) % THIS IS THE USEFUL METHOD (input images, variance and PSF to find streaks)
            
            import util.text.cs;
            
            I = []; % temporary storage for input image.
            
            % check if there are any images in the input varargin:
            if ~isempty(varargin) && isnumeric(varargin{1}) 
                I = varargin{1}; 
                if isscalar(varargin)
                    varargin = {};
                else
                    varargin = varargin(2:end);
                end
            else
                for ii = 1:2:length(varargin)
                    if cs(varargin{ii}, 'image')
                        I = varargin{ii+1};
                    end
                end
            end
            
            if isempty(I)
                return; % silently quit if there are no inputs.
            end
            
            obj.clear;
            obj.image = I;
            obj.parse(varargin{:}); % parse the rest of the inputs
            
            obj.run; % do the actual calculations. 
            
            if nargout>0
                streaks = obj.streak;
            end
            
        end
            
        function parse(obj, varargin) % recognize some keyword-value pairs as inputs and parameters. Does not recognize "image" as input! 
            
            import util.text.cs;
            import util.text.parse_bool;
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'variance')
                    obj.input_var = val;
                elseif cs(key, 'psf')
                    obj.input_psf = val;
                elseif cs(key, 'mean', 'use_mean')
                    obj.use_subtract_mean = parse_bool(val);
                elseif cs(key, 'mean', 'use_median')
                    obj.use_subtract_median = parse_bool(val);
                elseif cs(key, 'convolution', 'use_convolution')
                    obj.use_conv = parse_bool(val);
                elseif cs(key, 'crop', 'use_crop')
                    obj.use_crop = parse_bool(val);
                elseif cs(key, 'crop_size', 'size')
                    obj.crop_size = val;
                elseif cs(key, 'short', 'use_short')
                    obj.use_short = parse_bool(val);                    
                elseif cs(key, 'min_length', 'minimal_length', 'length')
                    obj.min_length = val;
                elseif cs(key, 'threshold')
                    obj.threshold = val;
                elseif cs(key, 'exclude', 'use_exclude')
                    obj.use_exclude = parse_bool(val);
                elseif cs(key, 'dx')
                    obj.exclude_dx = val;
                elseif cs(key, 'dy')
                    obj.exclude_dy = val;
                elseif cs(key, 'only_one', 'use_only_one')
                    obj.use_only_one = parse_bool(val);
                elseif cs(key, 'mean', 'use_mean')
                    obj.use_subtract_mean = parse_bool(val);
                elseif cs(key, 'filename')
                    obj.filename = val;
                elseif cs(key, 'batch_number')
                    obj.batch_num = val;
                elseif cs(key, 'frame_number')
                    obj.frame_num = val;
                elseif cs(key, 'section_number')
                    obj.section_num = val;
                elseif cs(key, 'offset_section', 'offset_original', 'corner')
                    obj.offset_original = val;
                elseif cs(key, 'was_convolved')
                    obj.was_convolved = parse_bool(val);
                elseif cs(key, 'original_image')
                    obj.original_image = val;
                end
                
            end
            
        end
        
        function run(obj)
            
            if obj.debug_bit
                fprintf('running streak detection on b= % 2d | f= % 3d | s= %d... ', obj.batch_num, obj.frame_num, obj.section_num);
            end

            if obj.use_preprocess
                obj.preprocess; % subtract background, remove point sources, get rid of NaN values
            else
                obj.image(isnan(obj.image)) = 0; % even when not preprocessing anything, must remove the NaN before going to FRT
            end
            
            obj.sendToFRT; % send to FRT, transpose, remove streaks, find best SNR
            
            if obj.debug_bit
                fprintf(' | length= %d | bestSNR= %f\n', round(obj.bestLength), obj.bestSNR);
            end
            
            if ~isempty(obj.streak)
                if obj.debug_bit, disp(['saving streaks with SNR= ' util.text.print_vec([obj.streak.snr])]); end
                obj.prev_streaks = [obj.prev_streaks obj.streak];
            end
            
            if obj.show_bit && obj.gui.check
                obj.show;
            end
            
        end
        
        function preprocess(obj)

            import util.img.crop2size;
            import util.stat.mean2;
            import util.stat.median2;
            
            if obj.use_crop % use this to crop down to power of 2
                [obj.image, gap] = crop2size(obj.image, obj.crop_size);
                obj.var_map = crop2size(obj.var_map, obj.crop_size);
                obj.offset_original = obj.offset_original + gap;
            end
            
            obj.im_size = util.vec.imsize(obj.image); 
            
            if obj.use_subtract_median
                obj.image = obj.image - median2(obj.image);
            end
            
            if obj.use_subtract_mean && obj.use_subtract_median==0
                obj.image = obj.image - mean2(obj.image);
            end

            obj.image(isnan(obj.image)) = 0; % If you cut the stars out and the were replaced with NaNs: subtract mean and only then replace NaNs with zeros.

            if obj.use_conv && ~isempty(obj.psf)
                obj.image = filter2(obj.psf, obj.image);
            end
            
        end
        
        function sendToFRT(obj)

            import util.stat.sum2;
            
            if 0 % I don't think we need any of this...
            % need these to normalize the final radon image
            obj.radon_var = obj.getRadonVariance(0); % if radon_var_map exists, use it. otherwise, get a radon_var_uni instead
            obj.radon_var_trans = obj.getRadonVariance(1);
            
            % move back to xy plane
            obj.radon_var = permute(obj.radon_var, [1,3,2]);
            obj.radon_var_trans = permute(obj.radon_var_trans, [1,3,2]);
            
            % calculate the geometric correction to the SNR
            N = pow2(nextpow2(obj.im_size(1)));
            th = atand((-N+1:N-1)./obj.im_size(1)); % angle range, in degrees
            G = max(abs(cosd(th)), abs(sind(th))); % geometric factor correction to the S/N
            NT = pow2(nextpow2(obj.im_size(2)));
            thT = atand((-NT+1:NT-1)./obj.im_size(2)); % angle range, in degrees
            GT = max(abs(cosd(thT)), abs(sind(thT))); % geometric factor correction to the S/N for the transposed image
            
            end
            
            R = radon.frt(obj.image, 'finder', obj, 'expand', obj.use_expand);
            temp_streak = obj.streak; % keep this here while looking for another streak in the transposed image
            obj.streak = radon.Streak.empty;
            
            RT = radon.frt(obj.image, 'finder', obj, 'expand', obj.use_expand, 'trans', 1);
            
            obj.streak = [temp_streak obj.streak]; % combine the two streaks
            
            if obj.use_only_one % if you get two streaks (one for each transpose) but only want the best one
                obj.streak = obj.bestStreak;
            end

            if 0 % I don't think we need any of this...
            obj.radon_image = R./sqrt(obj.radon_var.*obj.getNormFactorPSF.*G); % this normalization makes sure the result is in units of S/N, regardless of the normalization of the PSF. 
            obj.radon_image_trans = RT./sqrt(obj.radon_var_trans.*obj.getNormFactorPSF.*GT); % one makes sure the PSF is unity normalized, the other divides by the amount of noise "under" the PSF
            end
            
            for jj = 1:length(obj.streak)
                obj.streak(jj).was_convolved = obj.was_convolved || (obj.use_preprocess && obj.use_conv);                
            end
            
            % check if we want to save a copy of the subframe and input/final image in each of the streaks... 
            if obj.use_save_images
                for jj = 1:length(obj.streak)

%                     obj.streak(jj).image = obj.image; % we don't need this, it happens in "finalizeFRT"

%                     if obj.streak(jj).transposed % we don't need this, it happens in "finalizeFRT" ???
%                         obj.streak(jj).radon_image = obj.radon_image_trans;
%                     else
%                         obj.streak(jj).radon_image = obj.radon_image;
%                     end
                    
                    if ~isempty(obj.original_image)
                        obj.streak(jj).original_image = obj.original_image;
                        obj.streak(jj).offset_original = obj.offset_original;
                    end

                end
            else
                for jj = 1:length(obj.streak)
                    obj.streak(jj).image = [];
                    obj.streak(jj).subframe = [];
                    obj.streak(jj).radon_image = [];
                end
            end

%             new_snr = max(obj.last_snr, obj.bestSNR);

            obj.snr_values = [obj.snr_values obj.bestSNR];
            
        end
        
        function scan(obj, M, transpose) % do this on the final Radon image and on partial Radon matrices (for short streaks)
        % scans the given subsection of the FRT for streaks. 
        % M is the FRT matrix (can be full-FRT of subsection)
        % m is how many binary addition levels were performed to get it,
        % where we calculate m from the number of columns we have. 
            
            import util.stat.sum2;            
            import util.stat.maxnd;
            
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
            P = obj.getNormFactorPSF;
            
%             S = (size(M,3)+1)/2;
%             th = atand((-S+1:S-1)./S);
%             G = max(abs(cosd(th)), abs(sind(th))); % geometric factor correction to the S/N
            G = obj.getGeometricFactor(m);
            G = util.vec.topages(G);
            
            SNR = double(M)./sqrt(V.*P.*G); % this should be normalized to SNR units... 
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
            
            if isempty(obj.last_snr) || mx>obj.last_snr % this saves the best SNR regardless of transposition
                obj.last_snr = mx;
            end
            
            if mx>obj.threshold && (isempty(obj.streak) || mx>obj.streak.snr) % if a shorter streak is already saved (for this transpose), only update it if the new streak has higher SNR
                obj.streak = radon.Streak(obj, SNR, m, transpose, M(idx(1), idx(2), idx(3)), idx);
            end
            
        end
        
        function finalizeFRT(obj, radon_image, transpose) % call this at the end of the FRT function, if there is a finder

            m = log2(size(radon_image,2)+1)-1; % logarithmic step we are in
            V = obj.getRadonVariance(transpose, m);
            V = permute(V, [1,3,2]);
            P = obj.getNormFactorPSF;
            G = obj.getGeometricFactor(m);
            
            if transpose
                obj.radon_image_trans = radon_image./sqrt(V.*P.*G);
            else
                obj.radon_image = radon_image./sqrt(V.*P.*G);
            end
            
            if ~isempty(obj.streak) && obj.streak.radon_dx < obj.min_length % kill streaks that are too short (e.g. bright point-sources)
                obj.streak = radon.Streak.empty;
            end

            if ~isempty(obj.streak) % if we found a streak (that survived the min_length limit
                
                obj.streak.image = obj.image;
                obj.streak.radon_image = radon_image;
                
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
                obj.radon_var_uni = {};
                obj.radon_var_map = {};
            end
            
            if isempty(obj.input_var) || isscalar(obj.input_var)
                V = obj.radon_var_uni{transpose+1}{step}*obj.noise_var; % lazy load through getter function
            else
                V = obj.radon_var_map{transpose+1}{step}; % lazy load through getter function
            end
            
        end
                
        function purgeStreaks(obj) % depricated method (use as reference for sigma clipping and fit to EV distribution)
            
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
                if isempty(obj.streak)
                    return;
                end
                s = obj.streak(obj.getStreakIndex);
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
                N = length(obj.streak);
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
                    [~,val] = max([obj.streak.snr]);
                elseif obj.display_index>length(obj.streak) || obj.display_index<1
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
                    obj.display_index = length(obj.streak);
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
                N = length(obj.streak);
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

            obj.images = [];
            obj.radon_image = [];
            obj.radon_image_trans = [];
            obj.radon_var_uni = {};
            obj.radon_var_map = {};

            for ii = 1:length(obj.streak)
                obj.streak(ii).radon_image = [];
                obj.streak(ii).subframe = [];
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
            
            s_vec = [obj.streak obj.prev_streaks]; % go over all streaks
            
            for ii = 1:length(s_vec)
                
                s = s_vec(ii);
                
                I = s.image;
                
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

