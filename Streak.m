classdef Streak < handle
% This object describes a single streak that is found using the FRT, 
% including the raw results of the finding algorithm, and the streak 
% parameters derived from the results. 
% 
% Usually a radon.Finder object will have a vector of these objects saved 
% after it has finished going over some images. 
%
% The Streak object can optionally keep a copy of the original image and 
% the Radon image where it was detected, as well as housekeeping data such
% as the filename, batch number and frame number of that image. 
%
% The Streak object also has some functionality, such as plotting the image
% and Radon image with the streak parameters printed over it, as well as 
% producing a subtracted image (the image minus the streak itself). 
% 
    
    properties % outputs
        
        im_size; % size of original image (after transposed, if used. can be different than size(input_image))
        
        input_image; % subtracted image given to finder
        radon_image; % full Radon image for the correct transposed
        subframe; % subframe where the streak was located (for short streaks. Otherwise, the same as "radon_image")
        psf; % the PSF, after normalization, used for convolving the image (if any)
        
        % these help find the correct image where the streak exists
        frame_num = 1; % which frame in the batch
        batch_num = 1; % which batch in the run
        filename = ''; % which file it came from 
        
        threshold = []; % which threshold was used (if any)
        is_short = 0; % if yes, the subframe size is smaller than the full radon size
        
        % how bright the treak was
        I; % intensity (brightness per unit length)
        snr; % signal to noise ratio for this detection
        snr_fwhm;
        count; % number of photons in the whole frame
        
        % streak coordinates
        L; % length of streak (pixels)
        th; % angle in degrees (th=0 is on the x axis, positive angles is going down...)
        a; % y=ax+b
        b; % y=ax+b 
        x0; % point of crossing of streak (or continuation of streak) with x axis
        x1; % start point in x
        x2; % end point in x
        y1; % start point in y
        y2; % end point in y
        
        % radon coordinates
        transposed = 0; % did we transposed the image before doing the FRT (i.e. larger than 45 degrees)
        radon_step; % what logarithmic step "m" in the algorithm was this detection
        radon_max_idx; % index where in the Radon-subframe the maximum was located
        radon_y; % position coordinate in the Radon image (may be a transposed input image...)
        radon_dy; % slope coordinate in the Radon image (may be a transposed input image...)
        radon_y_var; % variance of the position
        radon_dy_var; % variance of the slope        
        radon_ydy_cov; % cross correlation of the slope-position errors
        radon_x1; % start of subsection of the original image (may be a transposedd input image)
        radon_x2; % end of subsection of the original image (may be a transposedd input image)
        
    end
    
    properties(Dependent=true)

        subframe_size; % size of subframe in Radon space (could be full radon size)
        radon_dx; % just x2-x1+1
        midpoint_x; % just (x1+x2)/2
        midpoint_y; % just (y1+y2)/2
        
    end
        
    properties % switches/controls
        
        noise_var = 1; % estimate of average noise variance 
        psf_sigma = 2; % estimate of PSF width parameter (assumed Gaussian)
        
        subtract_psf_widths = 3; % if not given by Finder, this tells "subtractStreak" how many PSF widths to remove around the streak position
        
        debug_bit = 1;
        
    end
        
    properties(Hidden=true)
       
        was_expanded = 0; % check if original image was expanded before FRT
        was_convolved = 1; % check if original image was convolved with PSF
        num_psfs_peak_region = 5; % rough estimate of the region around the Radon peak we want to cut for error estimates
        num_snr_peak_region = 2; % how many S/N units below maximum is still inside the peak region
        peak_region; % a map of the peak region, with only the part above the cut (peak-num_snr_peak_region) not zeroed out (used for error estimates)
        
        version = 1.03;
        
    end
    
    methods % constructor
        
        function obj = Streak(finder, subframe, log_step, transposed, count, index) % must give a bunch of parameters to create a legal streak (should be done only by radon.Finder)
            
            if nargin>0 && ~isempty(finder) && isa(finder, 'radon.Streak') % copy constructor
                
                if obj.debug_bit>1, fprintf('Streak copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(finder);
                
            else
                
                if obj.debug_bit>1, fprintf('Streak constructor v%4.2f\n', obj.version); end
                
                if nargin>=6 && ~isempty(finder)

                    if ~isa(finder, 'radon.Finder')
                        error('must input a radon.Finder object to constructor of radon.Streak');
                    end

                    if nargin<2 || isempty(subframe), error('must give a "subframe" to Streak constructor'); end
                    if nargin<3 || isempty(log_step), error('must specify the "log_step" for Streak constructor'); end
                    if nargin<4 || isempty(transposed), error('must specify the "transposed" for Streak constructor'); end
                    if nargin<5 || isempty(count), error('must specify the "count" for Streak constructor'); end
                    if nargin<6 || isempty(index), error('must specify the "index" for Streak constructor'); end
                                        
                    obj.subframe = subframe;
                    obj.radon_step = log_step;
                    obj.transposed = transposed;
                    obj.count = count;
                    obj.radon_max_idx = index;
                    obj.snr = subframe(index(1), index(2), index(3));
                    
                    try
                        obj.update_from_finder(finder); % get all the generic details straight from the finder object
                        obj.calculate; % go from raw Radon coordinates to streak parameters
                    catch ME
                        warning(ME.getReport);
                    end

                end
                
            end
            
        end
                
    end
    
    methods % reset/clear
        % these streaks are not reuseable
    end
    
    methods % getters
        
        function val = get.subframe_size(obj)
            
            val = size(obj.subframe);
            
        end
        
        function val = get.radon_dx(obj)
            
            val = obj.radon_x2 - obj.radon_x1 + 1;
            
        end
        
        function val = get.midpoint_x(obj)
           
            val = (obj.x1+obj.x2)/2;
            
        end
        
        function val = get.midpoint_y(obj)
           
            val = (obj.y1+obj.y2)/2;
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function update_from_finder(obj, finder) % get any parameters you can from the Finder object
            
            prop_list = properties(obj); 
            
            for ii = 1:length(prop_list) % load all properties that exist in "obj" and also in "finder"
                name = prop_list{ii};
                if isprop(finder, name)
                    obj.(name) = finder.(name); % includes: radon_image, psf, frame_num, batch_num, filename, threshold, noise_var, psf_sigma, debug_bit
                end
            end
            
            obj.was_convolved = finder.use_conv;
            obj.was_expanded = finder.use_expand;
            obj.im_size = finder.im_size_tr; % the size of the image that was actually entered into the FRT
                        
        end
        
        function calculate(obj, psf_sigma, noise_var) % transform raw Radon coordinates into streak parameters
            
            import util.stat.max2;
            import util.stat.sum2;
            
            if nargin>=2 && ~isempty(psf_sigma)
                obj.psf_sigma = psf_sigma;
            end
            
            if nargin>=3 && ~isempty(noise_var)
                obj.noise_var = noise_var;
            end
            
            obj.is_short = ndims(obj.subframe)>2 && size(obj.subframe,2)>1;
            obj.radon_x1 = (2.^obj.radon_step)*(obj.radon_max_idx(2)-1)+1;
            obj.radon_x2 = (2.^obj.radon_step)*(obj.radon_max_idx(2));
            obj.radon_y = obj.radon_max_idx(1);
            obj.radon_dy = obj.radon_max_idx(3) - ceil(size(obj.subframe,3)/2);
            
            % calculate everything assuming there was no transform...
            obj.x1 = obj.radon_x1;
            obj.x2 = obj.radon_x2;
            obj.a = obj.radon_dy/obj.radon_dx;
            offset = floor((obj.subframe_size(1)-obj.im_size(1))/2); % this difference is from using expandMatrix
            obj.b = obj.radon_y - offset - obj.radon_x1*obj.a;
            obj.x0 = -obj.b./obj.a;
            obj.th = atand(obj.a); % this theta is inverted when transposed==1
            obj.L = obj.radon_dx./cosd(obj.th);
            f = abs(cosd(obj.th));
            obj.I = obj.snr*sqrt(2.*sqrt(pi).*obj.psf_sigma.*obj.noise_var./(obj.L.*f));
            obj.snr_fwhm = obj.I*0.81./sqrt(obj.noise_var);
            obj.y1 = obj.a*obj.x1 + obj.b;
            obj.y2 = obj.a*obj.x2 + obj.b;
            
            % get the errors on the parameters...
            R = permute(obj.subframe(:, obj.radon_max_idx(2), :), [1,3,2]);
            xmax = obj.radon_max_idx(3);
            ymax = obj.radon_max_idx(1);
            
            S = round(obj.psf_sigma*obj.num_psfs_peak_region); % assume we only want a small section around the maximum
            x1 = xmax-S; if x1<1, x1 = 1; end
            x2 = xmax+S; if x2>size(R,2), x2 = size(R,2); end
            y1 = ymax-S; if y1<1, y1 = 1; end
            y2 = ymax+S; if y2>size(R,1), y2 = size(R,1); end
            C = R(y1:y2, x1:x2); % section around the max
            
            [X,Y] = meshgrid(1:size(C,2), 1:size(C,1));
            [mx,idx] = max2(C);
            
            C(C<mx-obj.num_snr_peak_region) = 0;
            
            X = X - idx(2);
            Y = Y - idx(1);
            
            obj.radon_dy_var = sum2(C.*X.^2)./sum2(C);
            obj.radon_y_var = sum2(C.*Y.^2)./sum2(C);
            obj.radon_ydy_cov = sum2(C.*X.*Y)./sum2(C);
            
            obj.peak_region = C;
            
            if obj.transposed % now what happens if the image was transformed?
                
                temp1 = obj.y1;
                temp2 = obj.y2;
                obj.y1 = obj.x1;
                obj.y2 = obj.x2;
                obj.x1 = temp1;
                obj.x2 = temp2;
                
                obj.a = 1./obj.a;
                obj.b = obj.x0;
                obj.x0 = -obj.b/obj.a;
                obj.th = 90 - obj.th;
                
            end
            
            % pixel values are rounded to nearest pixel
            obj.b = round(obj.b);
            obj.x1 = round(obj.x1);
            obj.x2 = round(obj.x2);
            obj.y1 = round(obj.y1);
            obj.y2 = round(obj.y2);
            obj.x0 = round(obj.x0);            
            
        end
        
        function M_sub = subtractStreak(obj, M_in, width) % subtract this streak (+error ellipse) from image
            
            if nargin<3 || isempty(width)
                width = obj.subtract_psf_widths;
            end
            
            % these are really rough estimates. Can improve this by looking at the error ellipse and subtracting y and dy values inside that range only 
            shift = -obj.psf_sigma*width:obj.psf_sigma*width;
%             dy = -obj.psf_sigma*width:obj.psf_sigma*width;
            
            M_sub = M_in;
            
            for ii = 1:length(shift)

                if obj.transposed
                    x1 = obj.x1 + shift(ii);
                    x2 = obj.x2 + shift(ii);
                    y1 = obj.y1;
                    y2 = obj.y2;
                else
                    x1 = obj.x1;
                    x2 = obj.x2;
                    y1 = obj.y1 + shift(ii);
                    y2 = obj.y2 + shift(ii);
                end
                
                [xlist, ylist] = radon.listPixels(x1,x2,y1,y2,size(M_in));

                idx = sub2ind(size(M_in), ylist{1}, xlist{1});
                M_sub(idx) = 0;

            end
            
        end
        
    end
    
    methods % plotting tools / GUI
          
        function val = getSubframeSlice(obj) % gives the radon-subframe where maximum happened, permuted to show image
           
            val = permute(obj.subframe(:,obj.radon_max_idx(2),:), [1,3,2]);
            
        end
        
        function val = getOriginalBand(obj) % shows the region of the input image in the corresponding radon-subframe of the streak
            
            if obj.transposed
                val = obj.input_image(obj.radon_x1:obj.radon_x2,:);
            else
                val = obj.input_image(:,obj.radon_x1:obj.radon_x2);
            end
            
        end
        
        function xywh = getBoxBoundaries(obj, width, height) % rectangle position vector around the streak in the original image
            
            if nargin<2 || isempty(width)
                width = obj.L.*2;
            end
            
            if nargin<3 || isempty(height)
                height = width;
            end
            
            x = round(obj.midpoint_x - width/2);
            y = round(obj.midpoint_y - height/2);
            w = round(width);
            h = round(height);
            
            if x<1, w = w-(1-x); x=1; end
            if y<1, h = h-(1-y); y=1; end
            
            xywh = [x y w h];            
            
        end
        
        function I = getCutout(obj, width) % a cutout around the streak in the original image
           
            if nargin<2 || isempty(width)
                width = obj.L.*2;
            end
                        
            xywh = obj.getBoxBoundaries(width, width);
            
            x1 = xywh(1);
            x2 = x1+xywh(3)-1;            
            y1 = xywh(2);
            y2 = y1+xywh(4)-1;
            
            if x1<1, x1 = 1; end
            if y1<1, y1 = 1; end            
            if x2>size(obj.input_image,2), x2 = size(obj.input_image,2); end            
            if y2>size(obj.input_image,1), y2 = size(obj.input_image,1); end
            
            I = util.img.pad2size(obj.input_image(y1:y2, x1:x2), width);
            
        end
        
        function show(obj, varargin) % plots both original and Radon image with streak parameters
            
            import util.plot.show;
            import util.plot.inner_title;
            import util.text.cs;
            import util.text.parse_bool;
            
            input_image = [];
            radon_image = [];
            font_size = [];
            use_publishable = [];
            parent = [];
                        
            % parse varargin
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'input_image')
                    input_image = val;
                elseif cs(key, 'radon_image')
                    radon_image = val;
                elseif cs(key, 'font_size')
                    font_size = val;
                elseif cs(key, {'use_publishable', 'publishable'})
                    use_publishable = parse_bool(val);
                elseif cs(key, 'parent')
                    parent = val;
                end
                
            end
                        
            if isempty(font_size)
                font_size = 24;
            end
            
            if isempty(input_image)
                input_image = obj.input_image;
            end
            
            if isempty(radon_image)
                radon_image = obj.radon_image;
            end
            
            if isempty(input_image) && isempty(radon_image)
                disp('cannot use "show" without original or final image!');
                return;
            end
            
            if isempty(parent)
                parent = gcf;
            end
                        
            if isvalid(parent)
                delete(parent.Children);
            end
                        
            ax1 = axes('Parent', parent, 'Position', [0.05 0.17 0.4 0.75]);
            
            if ~isempty(input_image) 
                
                obj.showOriginal(varargin{:}, 'axes', ax1, 'FontSize', font_size);
                
                if use_publishable
                    inner_title(ax1, '(a)', 'Position', 'corner');
                end
                
            end
            
            ax2 = axes('Parent', parent, 'Position', [0.5 0.1 0.4 0.8]);
            
            if ~isempty(radon_image)
                
                obj.showRadon(varargin{:}, 'axes', ax2);

                if use_publishable
                    inner_title(ax2, '(b)', 'Position', 'corner', 'FontSize', font_size);
                end
                
            end
            
            drawnow;
            
        end
           
        function showOriginal(obj, varargin) % show the original image, optionally with guidelines to show the found coordinates
            
            import util.plot.show;
            import util.plot.inner_title;
            import util.text.*;
            import util.stat.*;
            
            given_snr = [];
            given_x0 = [];
            given_y0 = [];
            given_th = [];
            given_I = [];

            line_offset = [];
            rect_size = [];
            font_size = [];    
            corner_title = '(a)';
            
            input_image = [];
           
            ax = [];
            
            use_monochrome = 0;
            use_publishable = 0;
            
            % parse varargin
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'snr')
                    given_snr = val;
                elseif cs(key, 'x0')
                    given_x0 = val;
                elseif cs(key, 'y0')
                    given_y0 = val;
                elseif cs(key, 'theta')
                    given_th = val;
                elseif cs(key, 'intensity')
                    given_I = val;
                elseif cs(key, 'streak')
                    streak = val;
                    if isa(streak, 'sim.Streak')
                        given_I = streak.I;
                        given_th = streak.th;
                        given_x0 = streak.x0;
                        given_snr = streak.getSNR;
                    end
                elseif cs(key, 'font_size')
                    font_size = val;
                elseif cs(key, 'line_offset')
                    line_offset = val;
                elseif cs(key, {'rectangle_size', 'rect_size'})
                    rect_size = val;
                elseif cs(key, 'input_image')
                    input_image = val;
                elseif cs(key, {'axis', 'axes'})
                    ax = val;
                elseif cs(key, {'monochrome', 'use_monochrome'})
                    use_monochrome = parse_bool(val);
                elseif cs(key, {'publishable', 'use_publishable'})
                    use_publishable = parse_bool(val);
                elseif cs(key, 'corner_title')
                    corner_title = val;
                end
                
            end
            
            if isempty(font_size)
                font_size = 24;
            end
            
            if use_publishable
                use_monochrome = 1;
            end
            
            if isempty(input_image)
                input_image = obj.input_image;
            end
            
            if use_monochrome                
                line_color = 'white';
            else
                line_color = 'red';
            end
                                    
            if isempty(ax)
                ax = gca;
            end
            
            show(input_image, 'axes', ax, 'autodyn', 1, 'monochrome', use_monochrome); 

            if use_publishable

                title(ax, '');

                inner_title(ax, corner_title, 'Position', 'corner', 'FontSize', font_size, 'FontColor', line_color);

            else

                if isempty(given_I)
                    title(ax, 'original image', 'FontSize', font_size);
                else
                    title(ax, ['original image (' f2s(given_I) 'e/pix)'], 'FontSize', font_size);
                end                    

                if ~use_publishable
                    obj.drawStats(ax, font_size);
                end

            end

            if ~isempty(line_offset)
                obj.drawGuidelines(ax, size(input_image), line_offset, line_color, use_publishable);
            end
            
            ax.FontSize = font_size;
            
            drawnow;
            
        end
        
        function drawStats(obj, ax, font_size) % plot some streak statistics on top of the image
            
            import util.plot.*;
            import util.text.*;
            
%             str = sprintf('b= %d | f= %d | x0= %d | \\theta= %3.0f^o | I= %3.1f', obj.batch_num, obj.frame_num, obj.x0, obj.th, obj.I);
            
            str = sprintf('b =%d | f= %d', obj.batch_num, obj.frame_num);
            if ~isempty(obj.batch_num) || ~isempty(obj.frame_num)
                inner_title(str, 'Position', 'Left', 'ax', ax, 'FontSize', font_size);
            end
            
            str = sprintf('x0= %d | \\theta = %3.0f^o', obj.x0, obj.th);
            inner_title(str, 'Position', 'Bottom', 'ax', ax, 'FontSize', font_size);
            
            str = sprintf('I= %3.1f | \\sigma= %3.1f | SNR= %3.1f', obj.I, sqrt(obj.noise_var), obj.snr);
            inner_title(str, 'Position', 'Right', 'ax', ax, 'FontSize', font_size);
            
        end
        
        function drawGuidelines(obj, ax, im_size, line_offset, line_color, use_publishable) % draw two dashed lines around the location of the streak in the original image
                        
                if line_offset==1
                    line_offset = 0.02*im_size(1);
                end

                if use_publishable
                    
                    N =3;
                    
                    line(ax, [obj.x1 obj.x1+line_offset*N*cosd(obj.th)]-line_offset, [obj.y1 obj.y1+line_offset*N*sind(obj.th)], 'Color', 'white', 'LineStyle','-', 'LineWidth', 2);
                    line(ax, [obj.x2-line_offset*N*cosd(obj.th) obj.x2]-line_offset, [obj.y2-line_offset*N*sind(obj.th) obj.y2], 'Color', 'white', 'LineStyle','-', 'LineWidth', 2);
                    
                    line(ax, [obj.x1 obj.x1+line_offset*N*cosd(obj.th)]+line_offset, [obj.y1 obj.y1+line_offset*N*sind(obj.th)], 'Color', 'white', 'LineStyle','-', 'LineWidth', 2);
                    line(ax, [obj.x2-line_offset*N*cosd(obj.th) obj.x2]+line_offset, [obj.y2-line_offset*N*sind(obj.th) obj.y2], 'Color', 'white', 'LineStyle','-', 'LineWidth', 2);
                    
                else
                    line(ax, [obj.x1 obj.x2]-line_offset, [obj.y1 obj.y2], 'Color', line_color, 'LineStyle','--', 'LineWidth', 2);
                    line(ax, [obj.x1 obj.x2]+line_offset, [obj.y1 obj.y2], 'Color', line_color, 'LineStyle','--', 'LineWidth', 2);
                end

        end
        
        function showRadon(obj, varargin) % show the Radon image 
            
            import util.plot.show;
            import util.plot.inner_title;
            import util.text.*;
            import util.stat.*;
            
            given_snr = [];
            given_x0 = [];
            given_y0 = [];
            given_th = [];
            given_I = [];

            rect_size = [];
            font_size = [];            
            expand = [];
            
            radon_image = [];
           
            ax = [];
            
            use_monochrome = 0;
            use_publishable = 0;
            
            % parse varargin
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'snr')
                    given_snr = val;
                elseif cs(key, 'x0')
                    given_x0 = val;
                elseif cs(key, 'y0')
                    given_y0 = val;
                elseif cs(key, 'theta')
                    given_th = val;
                elseif cs(key, 'intensity')
                    given_I = val;
                elseif cs(key, 'streak')
                    streak = val;
                    if isa(streak, 'sim.Streak')
                        given_I = streak.I;
                        given_th = streak.th;
                        given_x0 = streak.x0;
                        given_snr = streak.getSNR;
                    end
                elseif cs(key, 'font_size')
                    font_size = val;
                elseif cs(key, {'rectangle_size', 'rect_size'})
                    rect_size = val;
                elseif cs(key, 'radon_image')
                    radon_image = val;
                elseif cs(key, {'axis', 'axes'})
                    ax = val;
                elseif cs(key, {'monochrome', 'use_monochrome'})
                    use_monochrome = parse_bool(val);
                elseif cs(key, {'publishable', 'use_publishable'})
                    use_publishable = parse_bool(val);
                end
                
            end
            
            if isempty(expand)
                expand = 0;
            end
            
            if isempty(font_size)
                font_size = 24;
            end
            
            if use_publishable
                use_monochrome = 1;
            end
            
            if isempty(radon_image)
                radon_image = obj.radon_image;
            end
            
            if use_monochrome                
                line_color = 'white';
            else
                line_color = 'red';
            end
                                    
            if isempty(ax)
                ax = gca;
            end
            
            if expand && obj.im_size(1)==size(radon_image,1)
                radon_image = [NaN(size(radon_image)); radon_image; NaN(size(radon_image))];
            end
            
            size_y = obj.im_size(1);
            offset = size(radon_image,1)-2*size_y;
            max_dy = (size(radon_image,2)+1)/2;

            show(radon_image, 'axes', ax, 'xvalues', -max_dy:max_dy/4:max_dy, 'yvalues', -offset:100:size_y+offset, 'monochrome', use_monochrome);

            if obj.transposed
                inner_title(ax, 'x0', 'position', 'left', 'FontSize', font_size);
                inner_title(ax, '\Deltax', 'position', 'bottom', 'FontSize', font_size);
            else                  
                inner_title(ax, 'y0', 'position', 'left', 'FontSize', font_size);
                inner_title(ax, '\Deltay', 'position', 'bottom', 'FontSize', font_size);
            end
        
            if use_publishable
                inner_title(ax, '(b)', 'Position', 'corner', 'FontSize', font_size);
                title(ax, '');
            else

                if isempty(given_snr)
                    title(ax, ['FRT image | SNR= ' f2s(obj.snr)], 'FontSize', font_size);
                else
                    title(ax, ['FRT image | SNR= ' f2s(obj.snr) '(expect: ' f2s(given_snr) ')'], 'FontSize', font_size);
                end

            end

            if rect_size

%                 x_val = obj.radon_dy - rect_size/2;
%                 y_val = obj.radon_y - rect_size/2;

                if obj.transposed
                    y_val = obj.x0 - rect_size/2;
                    x_val = obj.im_size(1)./obj.a - rect_size/2;
                else
                    y_val = obj.b - rect_size/2;
                    x_val = obj.im_size(2).*obj.a - rect_size/2;
                end

                rectangle(ax, 'Position', [x_val y_val rect_size rect_size], 'EdgeColor', line_color);
                rectangle(ax, 'Position', [x_val-1 y_val-1 rect_size+2 rect_size+2]);

            end
            
            ax.FontSize = font_size;
            
            drawnow;
            
        end
        
    end    
    
    methods % utilities
        
        function saveas(obj_vec, filename, varname) % just save a streak (or vector of streaks) into a MAT file
                    
            ext  = '';
            if nargin<3 || isempty(varname)
                [~, varname, ext] = fileparts(filename);
            end
            
            if isempty(ext)
                filename = [filename '.mat'];
            end
            
            if obj_vec(1).debug_bit
                disp(['saving Streak objects "' varname '" in file: ' filename]);
            end
            
            save_struct.(varname) = obj_vec;
            
            save(filename, '-struct', 'save_struct', '-v7.3');
            
        end
        
    end
    
end

