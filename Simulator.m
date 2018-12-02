classdef Simulator < handle
% This object generates streaks with known position and intensity to use as
% tests for the radon.Finder. 
% It can generate defined streaks or random streaks. 
% It contains a radon.Finder object and by default feeds the images it 
% creates to the Finder and displays the resulting image with guidelines 
% and streak parameters (using the GUI). 
% 
% There are four ways to tell the Simulator what streak positions we want:
% (1) Using the "pickPoints" command: Click on points in the axes and those
%     are translated to lines. Whenever you right click or hit enter the 
%     line picker starts a new shape (each shape is made of at least two 
%     points. You can draw multiple shapes (each shape can be one line). 
%     To exit the picker, just press enter or right click twice (e.g. make 
%     an empty shape). 
% (2) Enter the (x,y) coordinates directly: use x1,x2,y1,y2 to define the 
%     streak start and end points, in normalized units (i.e. from 0 to 1). 
% (3) Enter "line" parameters: use "line_length", "line_theta", "line_x0", 
%     "line_a", "line_b", "line_midpoint_x" and "line_midpoint_y". 
%     There should be no problem entering them in any order, but notice the 
%     "anchor" key word chooses what point of the streak stays put while 
%     the angle and length are changed (use "start", "end" or "middle"). 
% (4) Randomly generate parameters: set use_random=1 and choose the ranges
%     (or a scalar for constant value) for all parameters. 
%     specifically: length, angle, midpoint_x and _y, as well as intensity 
%     can be chosen to vary uniformly on any range. Do not choose values 
%     outside the reasonable bounds (e.g. midpoints should be in [0,1]). 
%     Use "num_ranfom_streaks" to make multiple streaks. 
%     
% NOTE: make sure to select the correct settings for the Finder, in order 
%       to successfully find the streaks (e.g. "use_short" and for multi-
%       streak "use_recursive"). 
%
% The Simulator keeps several images, for the line production process: 
% (a) "image_line" is a single-pixel width line in an empty image. 
% (b) "image_conv" is the same, only convolved with the PSF of the image. 
%     If "use_psf=0" this image is the same as (a). 
% (c) "image_final" is the same, only with noise added. We use simple, 
%     Gaussian white noise. Choose the variance using "bg_noise_var" and 
%     control source noise with "use_source_noise". 
% (d) "psf" is the actual PSF image used to convolve the line. 
%     By default we use a width=2 Gaussian, but other widths can be set 
%     using "psf_sigma", and by setting "use_external_psf=1" this is 
%     replaced by the user input PSF (must give one!). 
%
% NOTE 2: the image size is set by "im_size" (scalar or 2-element vector). 
% NOTE 3: the "intensity" parameter sets the counts/unit length where the 
%         unit length is the side of a pixel. Diagonal lines will have 
%         somewhat brighter pixels because the streak is longer inside the 
%         pixel. 
%
% A Graphic User Interface (GUI) for this class is also included in the 
% sub-package +gui. It is invoked using obj.makeGUI. 
    
    properties (Transient=true)
        
        gui@radon.gui.SimGUI; % a GUI object for this class

    end
    
    properties % objects
        
        finder@radon.Finder; % use this to find streaks in the images produced. 
        
    end
    
    properties % outputs
        
        image_line; % single pixel width line
        image_conv; % convolve the line with the PSF
        image_final; % add noise to the line
        
        psf; % assume Gaussian PSF (can also be given by user, must set use_external_psf=1)
        
        num_pixels = 0; % number of pixels in which there is a line.
        
    end
    
    properties % switches/controls
        
        im_size = 512; % what size image to produce

        bg_noise_var = 1; % counts^2 per pixel
        use_source_noise = 1; % shot noise 
        
        use_psf = 1; % if you want to convolve the image with a PSF        
        psf_sigma = 2; % width prameter of the Gaussian PSF
        use_external_psf = 0; % if you don't want to generate a simple, symmetric Gaussian, must input a "psf" yourself
        
        use_random = 0; % generate random streak coordinates each time 
        num_random_streaks = 1; % how many random streaks to be made in each image
        intensity_range = [2 10]; % uniformly distributed
        length_range = [0.1 1]; % fractions of image size for radnom generation
        midpoint_x_range = [0.2 0.8]; % fractions of image size
        midpoint_y_range = [0.2 0.8]; % fractions of image size
        angle_range = [0 180]; % degrees
        
        use_finder = 1; % if you only want to make streak images, without directly sending them to the Finder
        
        use_update_finder = 1; % give correct parameters to finder. should most likely leave this on
        
        anchor = 'start'; % what to keep steady when changing "line_length" or "line_theta" or "line_midpoint_x" or "_y". Choose "start" or "end" or "middle" 
        
        display_what = 'final'; % which image to show when GUI is on. Choose "final" or "conv" or "line" 
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        num_lines; % how many lines are actually produced in the image (including lines with zero intensity)
        
        intensity; % per unit length (length of pixel side)--> xi in the paper
        line_length; % in normalized (image size) units, from 0 to 1
        line_theta; % in degrees measured down from the x axis
        line_midpoint_x; % in normalized units
        line_midpoint_y; % in normalized units
        line_x0; % in normalized units
        line_a; % slope parameter: a=tan(theta)
        line_b; % intercept parameter, in normalized units
        
    end
    
    properties(Hidden=true)
       
        % these are filled with default values at construction
        default_im_size;
        default_intensity;
        default_bg_noise_var;
        default_use_source_noise;
        
        default_use_psf;
        default_psf_sigma;
        default_use_external_psf;
        
        default_use_random;
        default_length_range;
        default_angle_range;
        default_midpoint_x_range;
        default_midpoint_y_range;
        
        intensity_private = 10; % this allows us to increase the number of lines and just copy the last "intensity" value to 0
        
        % real coordinates of the streak start/end points (normalized units)
        x1 = 0.1;
        y1 = 0;
        x2 = 0.9;
        y2 = 1;
        
        anchor_list = {'start', 'end', 'middle'};        
        display_what_list = {'final', 'conv', 'line'};
        
        version = 1.04;
        
    end
    
    methods % constructor
        
        function obj = Simulator(other)
            
            if nargin>0 && ~isempty(other) && isa(other, 'radon.Simulator') % copy-constructor
                if obj.debug_bit, fprintf('Simulator copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(other);
            else
            
                if nargin>0 && ~isempty(other) && isa(other, 'radon.Finder') % use this to give an existing Finder
                    obj.finder = other; 
                else
                    obj.finder = radon.Finder;
                    obj.finder.use_exclude = 0;
                    obj.finder.use_only_one = 1;
                    obj.finder.use_short = 0;
                end

                if obj.debug_bit, fprintf('Simulator constructor v%4.2f\n', obj.version); end

                util.oop.save_defaults(obj); % put variable "X" inside "default_X"

                obj.initialize;

            end
            
        end
        
    end
    
    methods % reset/clear
        
        function initialize(obj)
            
            util.oop.load_defaults(obj); % put variable "default_X" back to "X"
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            % get ready for a new image
            obj.image_line = [];
            obj.image_conv = [];
            obj.image_final = [];
            obj.psf = [];
            obj.num_pixels = 0;
            
        end
        
    end
    
    methods % getters
        
        function val = get.im_size(obj) % make sure the im_size output as a 2-element vector
            
            val = util.vec.imsize(obj.im_size);
            
        end
          
        function val = get.num_lines(obj)
           
            val = max([length(obj.x1), length(obj.x2), length(obj.y1), length(obj.y2)]);
            
        end
        
        function val = get.intensity(obj) 
           
            val = obj.intensity_private;
            
        end
        
        function val = get.line_length(obj) % normalized units
           
            if isempty(obj.x1) || isempty(obj.y1) || isempty(obj.x2) || isempty(obj.y2)
                val = [];
            else
                val = sqrt((obj.x2-obj.x1).^2+(obj.y2-obj.y1).^2);
            end
            
        end
        
        function val = line_length_pixels(obj) % pixel units
           
            if isempty(obj.line_length) || isempty(obj.im_size)
                val = [];
            else
                dx = (obj.x2-obj.x1)*obj.im_size(2);
                dy = (obj.y2-obj.y1)*obj.im_size(1);
                val = sqrt(dx.^2+dy.^2);
            end
            
        end
        
        function val = get.line_theta(obj) % degrees
           
            if isempty(obj.x1) || isempty(obj.y1) || isempty(obj.x2) || isempty(obj.y2)
                val = [];
            else
                val = atan2d(obj.y2-obj.y1, obj.x2-obj.x1);
            end
            
        end
        
        function val = get.line_midpoint_x(obj) % normalized units
            
            if isempty(obj.x1) || isempty(obj.x2)
                val = [];
            else
                val = (obj.x2+obj.x1)/2;
            end
            
        end
        
        function val = line_midpoint_x_pixels(obj) % pixel units
           
            if isempty(obj.line_midpoint_x) || isempty(obj.im_size)
                val = [];
            else
                val = obj.line_midpoint_x.*obj.im_size(2);
            end
            
        end
        
        function val = get.line_midpoint_y(obj) % normalized units
            
            if isempty(obj.y1) || isempty(obj.y2)
                val = [];
            else
                val = (obj.y2+obj.y1)/2;
            end
            
        end
        
        function val = line_midpoint_y_pixels(obj) % pixel units
           
            if isempty(obj.line_midpoint_y) || isempty(obj.im_size)
                val = [];
            else
                val = obj.line_midpoint_y*obj.im_size(1);
            end
            
        end
        
        function val = get.line_a(obj)
           
            if isempty(obj.x1) || isempty(obj.y1) || isempty(obj.x2) || isempty(obj.y2)
                val = [];
            else
                val = (obj.y2-obj.y1)./(obj.x2-obj.x1);
            end
            
        end
        
        function val = get.line_b(obj) % normalized units
           
            if isempty(obj.line_a)
                val = [];
            else
                val = obj.y1 - obj.line_a.*obj.x1;
            end
            
        end
        
        function val = line_b_pixels(obj) % pixel units
           
            if isempty(obj.line_b) || isempty(obj.im_size)
                val = [];
            else
                val = obj.line_b.*obj.im_size(1);
            end
            
        end
        
        function val = get.line_x0(obj)  % normalized units
           
            if isempty(obj.line_a)
                val = [];
            else
                val = -obj.line_b./obj.line_a;
            end
            
        end        
        
        function val = line_x0_pixels(obj) % pixel units
           
            if isempty(obj.line_x0) || isempty(obj.im_size)
                val = [];
            else
                val = obj.line_x0.*obj.im_size(2);
            end
            
        end
        
        function val = has_line(obj) % make sure the (x,y) coordinates are not empty
            
            if isempty(obj.x1) || isempty(obj.x2) || isempty(obj.y1) || isempty(obj.y2)
                val = 0;
            else
                val = 1;
            end
            
        end
        
        function val = isVertical(obj, idx) % for lines closer to vertical (>45 degrees)
            
            if nargin<2 
                idx = [];
            end

            val = (obj.line_theta>=45 & obj.line_theta<=135) |...
                    (obj.line_theta<=-45 & obj.line_theta>=-135) | ...
                    (obj.line_theta>=225 & obj.line_theta<315);
                
            if ~isempty(idx)
                val = val(idx);
            end
            
        end
        
    end
        
    methods % setters
                
        function set.num_lines(obj, val) % makes sure all parameters (x, y, intensity, num_pixels) are vectors of the same length
           
            if val<1
                val = 1;
            end
            
            if isempty(obj.intensity)
                obj.intensity = obj.default_intensity;
            end
            
            if isempty(obj.x1)
                obj.x1 = 0;
            end
            
            if isempty(obj.y1)
                obj.y1 = 0;
            end
            
            if isempty(obj.x2)
                obj.x2 = 1;
            end
            
            if isempty(obj.y2)
                obj.y2 = 1;
            end
            
            if length(obj.intensity_private)<val
                obj.intensity_private = [obj.intensity_private, ones(1, val-length(obj.intensity_private))*obj.intensity_private(end)];
            elseif length(obj.intensity_private)>val
                obj.intensity_private = obj.intensity_private(1:val);
            end
            
            if length(obj.x1)<val
                obj.x1 = [obj.x1 ones(1, val-length(obj.x1))*obj.x1(end)];
            elseif length(obj.x1)>val
                obj.x1 = obj.x1(1:val);
            end
            
            if length(obj.x2)<val
                obj.x2 = [obj.x2 ones(1, val-length(obj.x2))*obj.x2(end)];
            elseif length(obj.x2)>val
                obj.x2 = obj.x2(1:val);
            end
            
            if length(obj.y1)<val
                obj.y1 = [obj.y1 ones(1, val-length(obj.y1))*obj.y1(end)];
            elseif length(obj.y1)>val
                obj.y1 = obj.y1(1:val);
            end
            
            if length(obj.y2)<val
                obj.y2 = [obj.y2 ones(1, val-length(obj.y2))*obj.y2(end)];
            elseif length(obj.y2)>val
                obj.y2 = obj.y2(1:val);
            end
            
            if length(obj.num_pixels)>val
                obj.num_pixels = obj.num_pixels(1:val);
            end
            
        end
        
        function set.intensity(obj, val)
            
            obj.num_lines = length(val);
            
            obj.intensity_private = val;
            
        end
        
        function set.line_length(obj, L)
            
            import util.text.*;
            
            obj.num_lines = length(L);

            theta = obj.line_theta; % keep the original angle
            
            if cs(obj.anchor, 'start')            
            
                obj.x2 = obj.x1 + cosd(theta).*L;
                obj.y2 = obj.y1 + sind(theta).*L;
                
            elseif cs(obj.anchor, 'end')
            
                obj.x1 = obj.x2 - cosd(theta).*L;
                obj.y1 = obj.y2 - sind(theta).*L;
                
            elseif cs(obj.anchor, 'middle')
                
                x = obj.line_midpoint_x;
                y = obj.line_midpoint_y;

                obj.x1 = x - cosd(theta).*L/2;
                obj.y1 = y - sind(theta).*L/2;
                obj.x2 = x + cosd(theta).*L/2;
                obj.y2 = y + sind(theta).*L/2;

            end
            
        end
        
        function set.line_theta(obj, theta)
            
            import util.text.*;
            
            obj.num_lines = length(theta);
            
            L = obj.line_length; % keep the original length
            
            if cs(obj.anchor, 'start')            
            
                obj.x2 = obj.x1 + cosd(theta).*L;
                obj.y2 = obj.y1 + sind(theta).*L;
                
            elseif cs(obj.anchor, 'end')
            
                obj.x1 = obj.x2 - cosd(theta).*L;
                obj.y1 = obj.y2 - sind(theta).*L;
                
            elseif cs(obj.anchor, 'middle')
                
                x = obj.line_midpoint_x;
                y = obj.line_midpoint_y;

                obj.x1 = x - cosd(theta).*L/2;
                obj.y1 = y - sind(theta).*L/2;
                obj.x2 = x + cosd(theta).*L/2;
                obj.y2 = y + sind(theta).*L/2;
                
            end
            
        end
        
        function set.line_midpoint_x(obj, val)

            obj.num_lines = length(val);

            x_old = obj.line_midpoint_x;
            x1 = obj.x1 + val - x_old;
            x2 = obj.x2 + val - x_old;
            obj.x1 = x1;
            obj.x2 = x2;
            
        end
        
        function set.line_midpoint_y(obj, val)

            obj.num_lines = length(val);

            y_old = obj.line_midpoint_y;
            
            y1 = obj.y1 + val - y_old;
            y2 = obj.y2 + val - y_old;
            
            obj.y1 = y1;
            obj.y2 = y2;
            
        end
        
        function set.line_x0(obj, val)
            
            obj.num_lines = length(val);

            old_x0 = obj.line_x0;
            
            obj.x1 = obj.x1 + val - old_x0;
            obj.x2 = obj.x2 + val - old_x0;
            
        end
        
        function set.line_a(obj, val)
            
            obj.line_theta = atand(val);
            
        end
        
        function set.line_b(obj, val)
            
            obj.num_lines = length(val);
            
            old_b = obj.line_b;
            
            obj.y1 = obj.y1 + val - old_b;
            obj.y2 = obj.y2 + val - old_b;
                        
        end
        
        function set.x1(obj, val)
           
            obj.x1 = util.vec.torow(val);
            
        end
        
        function set.x2(obj, val)
           
            obj.x2 = util.vec.torow(val);
            
        end
        
        function set.y1(obj, val)
           
            obj.y1 = util.vec.torow(val);
            
        end
        
        function set.y2(obj, val)
           
            obj.y2 = util.vec.torow(val);
            
        end
        
    end
    
    methods % calculations
        
        function val = trig_factor(obj, idx) % either sin(th) or cos(th), whichever is closer to 1
            
            if nargin<2 || isempty(idx)
                idx = [];
            end
            
            if isempty(idx)            
                val = max(abs(cosd(obj.line_theta)), abs(sind(obj.line_theta)));
            else
                val = max(abs(cosd(obj.line_theta(idx))), abs(sind(obj.line_theta(idx))));
            end
            
        end
        
        function val = calcSNR(obj) % what is the theoretical S/N for the given streaks
            
            L = obj.num_pixels./obj.trig_factor;
            
            I = obj.intensity; % already includes geometric factor... 
            
%             val = I.*sqrt(L.*obj.trig_factor./obj.bg_noise_var);
            val = I.*sqrt(L./obj.bg_noise_var);
            
            if obj.use_psf
                val = val./sqrt(2*sqrt(pi)*obj.psf_sigma);
            end
            
            val = abs(val);
            
        end
                
        function randomize(obj) % generate random streaks from uniform distribution of parameters
           
            L = [];
            th = [];
            mid_x = [];
            mid_y = [];
            I = [];
            
            for ii = 1:obj.num_random_streaks
            
                if length(obj.length_range)>1
                    L(ii) = rand*(obj.length_range(2)-obj.length_range(1))+obj.length_range(1);
                elseif isscalar(obj.length_range)
                    L(ii) = obj.length_range;            
                end

                if length(obj.angle_range)>1
                    th(ii) = rand*(obj.angle_range(2)-obj.angle_range(1))+obj.angle_range(1);
                elseif isscalar(obj.angle_range)
                    th(ii) = obj.angle_range;
                end

                if length(obj.midpoint_x_range)>1
                    mid_x(ii) = rand*(obj.midpoint_x_range(2)-obj.midpoint_x_range(1))+obj.midpoint_x_range(1);
                elseif isscalar(obj.midpoint_x_range)
                    mid_x(ii) = obj.midpoint_x_range;
                end

                if length(obj.midpoint_y_range)>1
                    mid_y(ii) = rand*(obj.midpoint_y_range(2)-obj.midpoint_y_range(1))+obj.midpoint_y_range(1);
                elseif isscalar(obj.midpoint_y_range)
                    mid_y(ii) = obj.midpoint_y_range;
                end
                
                if length(obj.intensity_range)>1
                    I(ii) = rand*(obj.intensity_range(2)-obj.intensity_range(1))+obj.intensity_range(1);
                elseif isscalar(obj.intensity_range)
                    I(ii) = obj.intensity_range;
                end                
                
%                fprintf('x= %f | y= %f | L= %f | th= %f\n', obj.line_midpoint_x, obj.line_midpoint_y, obj.line_length, obj.line_theta);
            
            end
            
            obj.line_length = L;
            obj.line_theta = th;
            obj.line_midpoint_x = mid_x;
            obj.line_midpoint_y = mid_y;
            obj.intensity = I;
            
        end
        
        function update_finder(obj) % give the Finder object the real variance and PSF 
            
%             if obj.finder.noise_var~=obj.bg_noise_var
%                 obj.finder.makeVarMap(obj.bg_noise_var);
%             end
            
            obj.finder.input_var = obj.bg_noise_var;
            obj.finder.input_psf = obj.psf;
                        
        end
        
        function run(obj) % generate streak(s) and optionally pass them to the Finder
           
            if obj.use_random
                obj.randomize;
            end
            
            obj.clear;
            obj.makeImage;
            
            if obj.use_finder && ~isempty(obj.finder)
                
                if obj.use_update_finder
                    obj.update_finder;
                end
                
                obj.finder.input(obj.image_final);
                
            end
            
            if ~isempty(obj.gui) && obj.gui.check
                obj.show('margin', 0.05);   
            end
            
        end
        
        function makeImage(obj) % fills the three images in turn: line, conv, final
            
            obj.drawLine;
            obj.runConv;
            obj.addNoise;
            
        end
        
        function [x_list, y_list] = listPixels(obj) % lists of x,y pixel coordinates of the streak(s) (uses "radon.listPixels")
            
            if isempty(obj.im_size)
                error('cannot get list of pixels without image size!');
            end
            
            if isempty(obj.x1) || isempty(obj.x2) || isempty(obj.y1) || isempty(obj.y2)
                error('Cannot get list of pixels without defining a line start/end point...');
            end
            
            x1 = obj.x1.*obj.im_size(2);
            x2 = obj.x2.*obj.im_size(2);
            y1 = obj.y1.*obj.im_size(1);
            y2 = obj.y2.*obj.im_size(1);
            
            [x_list, y_list, obj.num_pixels] = radon.listPixels(x1,x2,y1,y2,obj.im_size); % the radon.listPixels is used in other places as well
            
        end
                
        function drawLine(obj) % takes the list of x,y pixel coordinates and puts "intensity" in each pixel of the line
            
            obj.image_line = zeros(obj.im_size);
                
            [x_list, y_list] = obj.listPixels;
            
            I = obj.intensity;
            
            if length(I)<length(x_list)
                I = [I ones(1,(length(x_list)-length(I)))*I(end)];
            end
            
            for ii = 1:length(x_list)
                
                idx = sub2ind(obj.im_size, y_list{ii}, x_list{ii});
                
                if isempty(idx)
                    continue;
                end
                
                obj.image_line(idx) = obj.image_line(idx) + I(ii)./obj.trig_factor(ii); % intensity is per unit length (diagonal lines have a longer length in each pixel...)
                
            end
            
        end
        
        function runConv(obj) % convolve the image with the PSF
            
            if obj.use_psf
                
                if ~obj.use_external_psf % automatically generate a PSF (deletes user-input PSF)
                    k = util.img.gaussian2(obj.psf_sigma);
                    obj.psf = k./util.stat.sum2(k);
                end
                
                obj.image_conv = filter2(obj.psf, obj.image_line); % don't think we need util.fft.conv_f, the PSF is too small for FFT convolution. 
                
            else
                obj.image_conv = obj.image_line; % if no PSF is used, just copy "image_line"
            end
            
        end
        
        function addNoise(obj) % add white Gaussian noise
            
            if obj.bg_noise_var>0
                obj.image_final = normrnd(zeros(obj.im_size), sqrt(obj.bg_noise_var));
            else
                obj.image_final = zeros(obj.im_size);
            end
            
            if obj.use_source_noise
                obj.image_final = obj.image_final + normrnd(obj.image_conv, sqrt(obj.image_conv));
            else
                obj.image_final = obj.image_final + obj.image_conv;
            end
                        
        end
        
    end
    
    methods % plotting tools / GUI
        
        function show(obj, varargin) % display the image with streak(s) on it, and some statistics
            
            import util.text.cs;
            import util.text.parse_bool;
            
            margin = []; % add margin to be able to draw lines outside the frame
            ax = []; % by default us the GUI axes if it exists, otherwise gca
            font_size = 24;
            use_publishable = 0;
            
            for ii = 1:2:length(varargin)
               
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, 'margin')
                    margin = val;
                elseif cs(key, {'axes', 'axis'})
                    ax = val;
                elseif cs(key, 'font_size')
                    font_size = val;
                elseif cs(key, 'publishable', 'use_publishable')
                    use_publishable = parse_bool(val);
                end
                
            end
            
            if isempty(margin)
                margin = 0;
            end
            
            if isempty(ax)
                if ~isempty(obj.gui) && obj.gui.check
                    ax = obj.gui.axes_image;
                else
                    ax = gca;
                end
            end
            
            if cs(obj.display_what, 'final')
                I = obj.image_final;
            elseif cs(obj.display_what, 'conv')
                I = obj.image_conv;
            elseif cs(obj.display_what, 'line')
                I = obj.image_line;
            end
            
            if margin>0
                margin_pix = ceil(margin*obj.im_size);
                I = util.img.pad2size(I, obj.im_size+margin_pix*2);
                xval = -margin_pix(2)+2:obj.im_size(2)+margin_pix(2)+1;
                yval = -margin_pix(1)+2:obj.im_size(1)+margin_pix(1)+1;
                util.plot.show(I, varargin{:}, 'xvalues', xval, 'yvalues', yval);
            else
                util.plot.show(I, varargin{:});
            end
            
            % get the Finder results and show them also
            if obj.use_finder && ~isempty(obj.finder) && ~isempty(obj.finder.streak)
                
                for ii = 1:length(obj.finder.streak)
                    obj.finder.streak(ii).drawGuidelines(ax, obj.im_size, 15, 'White', use_publishable);
                end
                
                obj.finder.streak(1).drawStats(ax, font_size);
                
            else
                delete(findobj(ax, 'type', 'line'));
            end
            
            if obj.has_line
                util.plot.inner_title(['SNR(calc)= ' util.text.f2s(obj.calcSNR)], 'ax', ax, 'Position', 'Top','FontSize', font_size-10);
            end
            
            if ~isempty(obj.gui)
                obj.gui.update;
            end
            
        end
        
        function makeGUI(obj) % create a GUI object (if it is empty) and make a GUI
            
            if isempty(obj.gui)
                obj.gui = radon.gui.SimGUI(obj);
            end
            
            obj.gui.makeGUI;
            
        end
        
        function pickPoints(obj, varargin) % select points on the axes to make them into lines
            
            import util.text.*;
            
            ax = [];
            
            for ii = 1:2:length(varargin)
                
                key = varargin{ii};
                val = varargin{ii+1};
                
                if cs(key, {'axis', 'axes'})
                   ax = val; 
                end
            
            end
            
            if isempty(ax)
                if ~isempty(obj.gui) && obj.gui.check
                    ax = obj.gui.axes_image;
                else
                    ax = gca;
                end
            end
            
            obj.show('margin', 0.05);

            x1 = [];
            x2 = [];
            y1 = [];
            y2 = [];
            
            for ii = 1:100
            
                [x,y] = getline(ax);

                if length(x)<2
                    break;
                end

                x = x./obj.im_size(2);
                y = y./obj.im_size(1);
                
                for jj = 1:length(x)-1
                    
                    x1(end+1) = x(jj);
                    x2(end+1) = x(jj+1);
                    y1(end+1) = y(jj);
                    y2(end+1) = y(jj+1);
                    
                end
            
            end
            
            obj.num_lines = length(x1);
            
            obj.x1 = x1;
            obj.x2 = x2;
            obj.y1 = y1;
            obj.y2 = y2;
            
            obj.makeImage;
            
            if ~isempty(obj.finder)
                
                if obj.use_update_finder
                    obj.update_finder;
                end
                
                obj.finder.input(obj.image_final);
                
            end
            
            obj.show('margin', 0.05);
            
        end
        
    end    
    
end
