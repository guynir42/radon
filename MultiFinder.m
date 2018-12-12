classdef MultiFinder < handle
% Scan an input image for streaks, using multiple thresholds and multiple 
% iterations on each image, so it can find multiple streaks and artefacts, 
% even in fairly bad subtraction images. 
%
% This object is meant to be used at a level above the Finder, it produces 
% a Finder object that works on each transposition and each threshold level. 
% In each iteration, streaks that are detected are removed from the image. 
% The threshold is logarithmically reduced, truncating the pixel values to 
% the new threshold, before running the Finder again. 
%
% The image can be sectioned into smaller pieces, so each piece is searched 
% iteratively over thresholds and streaks. 
%
    
    properties(Transient=true)
        
        gui@radon.gui.MultiGUI;
        
    end
    
    properties % objects
        
        finder@radon.Finder; 
%         prof@radon.Profiler;
        class@radon.Classifier;
        
        streaks@radon.Streak; % all streaks saved up to the last "reset"
        asteroids@radon.Streak; % everything that was classified as NEO
        satellites@radon.Streak; % everything that was classified as a satellite (should we split this to LEO objects?)
        artefacts@radon.Streak; % everything else that is not astronomical
        
        streaks_all@radon.Streak; % all streaks saved up to the last "reset"
        asteroids_all@radon.Streak; % everything that was classified as NEO
        satellites_all@radon.Streak; % everything that was classified as a satellite (should we split this to LEO objects?)
        artefacts_all@radon.Streak; % everything else that is not astronomical
        
    end
    
    properties % inputs/outputs
        
        image;
        image_subtracted; % full frame, non filtered image, with streaks removed
        image_processed; % current image or section that undergoes processing and streak subtractions
        image_output; % image after truncating bright pixels and removing streaks
        image_sections = {}; % cell used to store sections of the image
        offset_sections = {};
        
        input_var;
        var_sections = {};
        
    end
    
    properties % switches/controls
        
        threshold = 10;
        
        use_scan_thresh = 1;
        use_conv = 1;
        iter_max = 5; % can increase this I think
        
        use_sections = 0;
        section_size = 1024;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        input_psf;
        psf;
        psf_sigma;
        var_scalar;
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = MultiFinder(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'radon.MultiFinder')
                if obj.debug_bit, fprintf('MultiFinder copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('MultiFinder constructor v%4.2f\n', obj.version); end

                util.oop.save_defaults(obj); % put default values for "X" into "default_X"

                obj.finder = radon.Finder;
                obj.finder.use_conv = 0;
                obj.finder.use_short = 1;
                obj.finder.use_only_one = 0;
                
                obj.class = radon.Classifier;
                
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.finder.reset;
%             obj.prof.reset;
              
            obj.streaks_all = radon.Streak.empty;
            obj.asteroids_all = radon.Streak.empty;
            obj.satellites_all = radon.Streak.empty;
            obj.artefacts_all = radon.Streak.empty;

            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.image = [];
            obj.image_processed = [];
            obj.image_sections = {};
            obj.offset_sections = {};
            obj.var_sections = {};
            
            obj.streaks = radon.Streak.empty;
            obj.asteroids = radon.Streak.empty;
            obj.satellites = radon.Streak.empty;
            obj.artefacts = radon.Streak.empty;
            
        end
        
    end
    
    methods % getters
        
        function val = get.input_psf(obj)
            
            val = obj.finder.input_psf;
            
        end
        
        function val = get.psf(obj)
            
            val = obj.finder.psf;
            
        end
        
        function val = get.psf_sigma(obj)
            
            val = obj.finder.psf_sigma;
            
        end
        
        function val = get.var_scalar(obj)
            
            if isempty(obj.input_var)
                val = [];
            elseif isscalar(obj.input_var)
                val = obj.input_var;
            else
                val = util.stat.median2(obj.input_var);
            end
            
        end
        
    end
    
    methods % setters
        
        function set.input_psf(obj, val)
            
            obj.finder.input_psf = val;
            
        end
        
        function set.image(obj, val)
            
            obj.image = val;
            obj.image_subtracted = val;
            
        end
        
    end
    
    methods % calculations
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('image', []);
            input.input_var('variance', [], 'input_variance');
            input.input_var('psf', [], 'input_psf');
            input.input_var('filename', '');
            input.input_var('batch_num', []);
            input.input_var('frame_num', []);
            input.scan_vars(varargin{:});
            
            if isempty(input.image)
                return; % quietly quit if there are no images
            end
            
            obj.clear;
            
            obj.image = input.image;
            
            if ~isempty(input.psf)
                obj.input_psf = input.psf;
            end
            
            if ~isempty(input.variance)
                obj.input_var = input.variance;
            end
            
            if obj.use_sections
                
                obj.makeSections;
                
                for ii = 1:length(obj.image_sections)
                    
                    obj.image_processed = obj.preprocess(obj.image_sections{ii});
                    obj.runSection('offset', obj.offset_sections{ii}, 'section', ii,...
                        'variance', obj.var_sections{ii}, 'psf', input.psf, ...
                        'filename', input.filename, 'batch', input.batch_num, 'frame', input.frame_num);
                    
                end
                
            else
                
                obj.image_processed = obj.preprocess(obj.image);
                obj.runSection('variance', obj.input_var, 'psf', input.psf,...
                    'filename', input.filename, 'batch', input.batch_num, 'frame', input.frame_num);
                
            end
            
            % add final view of full image and streaks with types
            
%             obj.calcProfiles;
%             obj.calcClassifications;

            obj.streaks_all = [obj.streaks_all obj.streaks];
            obj.asteroids_all = [obj.asteroids_all obj.asteroids];                
            obj.satellites_all = [obj.satellites_all obj.satellites];                
            obj.artefacts_all = [obj.artefacts_all obj.artefacts];
            
        end
        
        function makeSections(obj, varargin)
            
            S = size(obj.image);
            
            axis1 = 1:obj.section_size:S(1);
            axis2 = 1:obj.section_size:S(2);
            
            obj.image_sections = {};
            obj.offset_sections = {};
            obj.var_sections = {};
            
            for ii = 1:length(axis1)
                for jj = 1:length(axis2)
                
                    if ii<length(axis1)
                        end_point1 = axis1(ii+1)-1;
                    else
                        end_point1 = S(1);
                    end
                    
                    if jj<length(axis2)
                        end_point2 = axis2(jj+1)-1;
                    else
                        end_point2 = S(2);
                    end
                
                    obj.image_sections{end+1} = obj.image(axis1(ii):end_point1, axis2(jj):end_point2);
                    obj.offset_sections{end+1} = [axis1(ii), axis2(jj)] - 1; % need to add this to y,x position in section to translate to full image
                    
                    if ~isscalar(obj.input_var)
                        obj.var_sections{end+1} = obj.input_var(axis1(ii).end_point1, axis2(jj):end_point2);
                    else
                        obj.var_sections{end+1} = obj.input_var;
                    end
                    
                end
            end
            
        end
        
        function I_out = preprocess(obj, I)
        
            I_out = I;
            
            if obj.use_conv
                I_out = filter2(obj.psf, I_out);
            end
            
        end
        
        function runSection(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('variance',[]);
            input.input_var('filename', '');
            input.input_var('batch_num', []);
            input.input_var('frame_num', []);
            input.input_var('section', 1, 'section_number');
            input.input_var('offset', [0 0]);
            input.scan_vars(varargin{:});
            
            I = obj.image_processed;
            V = input.variance; % scalar value of the variance (or median of the var map)
            
            if isempty(V)
                V = 1; 
            elseif ismatrix(V)
                V = util.stat.median2(V);
            end
                        
            if obj.use_scan_thresh
                log_values = log2(obj.threshold):log2(util.stat.max2(I)./sqrt(V));
                thresh_values = flip(2.^log_values);
            else
                thresh_values = obj.threshold;
            end
            
            for ii = 1:length(thresh_values)
                
                if obj.gui.check
                    obj.showProcessed('thresh', thresh_values(ii), 'section', input.section);
                    drawnow;
                end
                
                obj.finder.threshold = thresh_values(ii);
                
                mask = I./sqrt(V)>thresh_values(ii);
                
                I(mask) = thresh_values(ii).*sqrt(V);
                
                obj.image_output = I;
                
                if obj.debug_bit, fprintf('  Threshold= %4.2f | mask= % 6d (%8.6f%%)\n', thresh_values(ii), nnz(mask), nnz(mask)./numel(I)*100); end 
                
                temp_streaks = radon.Streak.empty;
                
                for jj = 1:obj.iter_max
                    
                    new_streaks = obj.finder.input('image', I, 'variance', input.variance,...
                         'was_convolved', obj.use_conv, 'original_image', obj.image, ...
                         'filename', input.filename, 'batch_num', input.batch_num,...
                         'frame_num', input.frame_num, 'section_num', input.section, 'offset', input.offset);
                    
                    if isempty(new_streaks)
                        break;
                    end
                    
                    for kk = 1:length(new_streaks)
                        
                        new_streaks(kk).prof.input(obj.image_subtracted);
                        obj.class.input(new_streaks(kk));
                        
                        I = new_streaks(kk).subtractStreak('image', I, 'replace', NaN);
                        obj.image_subtracted = new_streaks(kk).subtractStreak('image', obj.image_subtracted, 'replace', NaN, 'offset', new_streaks(kk).offset_original);
                        
                        % add the new streaks to the list
                        obj.streaks = [obj.streaks new_streaks(kk)];
                        
                        if new_streaks(kk).is_asteroid
                            obj.asteroids = [obj.asteroids new_streaks(kk)];
                        end
                        
                        if new_streaks(kk).is_satellite
                            obj.satellites = [obj.satellites new_streaks(kk)];
                        end
                        
                        if new_streaks(kk).is_artefact
                            obj.artefacts = [obj.artefacts new_streaks(kk)];
                        end
                        
                    end
                    
                    obj.image_output = I;
                    
                    % if we found any new streaks
                    temp_streaks = [temp_streaks new_streaks];
                    
                    if obj.gui.check
                        obj.showProcessed('thresh', thresh_values(ii),'section', input.section);
                        drawnow;
                    end
                    
                end
                
                if obj.debug_bit
                    if ~isempty(temp_streaks)
                        fprintf('  Found streaks with SNR: %s\n', util.text.print_vec([temp_streaks.snr])); 
                    end
                end
                
            end
            
        end
        
        function calcProfiles(obj)
            
            for ii = 1:length(obj.streaks)
                
                I = obj.image;
                
                for jj = 1:length(obj.streaks)
                    
                    if ii==jj, continue; end
                    
                    I = obj.streaks(jj).subtractStreak('image', I, 'replace', NaN, 'offset', obj.streaks(jj).offset_original);
                    
                end
                
                if obj.gui.check && 0
                    
                    util.plot.show(I, 'bias', 0, 'dyn', 10, 'ax', obj.gui.axes_image);
                    xlabel(obj.gui.axes_image, ['ii= ' num2str(ii)]);
                    obj.streaks(ii).drawGuidelines('ax', obj.gui.axes_image, 'size', size(I), 'offset', obj.streaks(ii).offset_original);
                    pause(2);
                    
                end
                
                obj.streaks(ii).prof.input('image', I);
                    
            end
            
        end
        
        function calcClassifications(obj)
            
            
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function makeGUI(obj) % fills the GUI object if needed
           
            if isempty(obj.gui)
                obj.gui = radon.gui.MultiGUI(obj);
            end
            
            obj.gui.make;
            
        end
        
        function show(obj, varargin)
            
            
            
        end
        
        function showOriginal(obj, varargin)
                        
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis');
            input.input_var('section', [], 'section_number');
            input.scan_vars(varargin{:});
            
            if isempty(input.axes)
                if obj.gui.check
                    input.axes = obj.gui.axes_image;
                else
                    input.axes = gca;
                end 
            end
            
            util.plot.show(obj.image, varargin{:});
            
            for ii = 1:length(obj.streaks)
                obj.streaks(ii).drawGuidelines('ax', input.axes, 'offset', obj.streaks(ii).offset_original);
            end
            
        end
        
        function showProcessed(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('axes', [], 'axis');
            input.input_var('threshold', []);
            input.input_var('section', 1, 'section_number');
            input.scan_vars(varargin{:});
            
            if isempty(input.axes)
                if obj.gui.check
                    input.axes = obj.gui.axes_image;
                else
                    input.axes = gca;
                end 
            end
            
            util.plot.show(obj.image_processed, 'axes', input.axes, 'bias', 0, 'dyn', input.threshold);
            
            num_sections = 1;
            if ~isempty(obj.image_sections) && obj.use_sections
                num_sections = length(obj.image_sections);
            end
            
            counter = 0;
            
            for ii = 1:length(obj.streaks)
                
                if input.section==obj.streaks(ii).section_num
                    
                    obj.streaks(ii).drawGuidelines('ax', input.axes, 'size', size(obj.image_processed), 'line_dist', 15);
                    counter = counter + 1;
                end
                
            end
            
            xlabel(input.axes, sprintf('section= %d/%d | thresh= %4.2f | streaks= %d', input.section, num_sections, input.threshold, counter));
            
        end
        
    end    
    
end

