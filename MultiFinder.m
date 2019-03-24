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
        
        second_fig;
        
    end
    
    properties % objects
        
        finder@radon.Finder; 
%         prof@radon.Profiler;
        class@radon.Classifier;
        
        streaks@radon.Streak; % all streaks saved up to the last "reset"
        asteroids@radon.Streak; % everything that was classified as NEO
        satellites@radon.Streak; % everything that was classified as a satellite (astronomical but too glinty)
        artefacts@radon.Streak; % everything else that is not astronomical
        too_longs@radon.Streak; % anything too long to get further analysis (image edges, artefacts, LEO satellites)
        
        streaks_all@radon.Streak; % all streaks saved up to the last "reset"
        asteroids_all@radon.Streak; % everything that was classified as NEO
        satellites_all@radon.Streak; % everything that was classified as a satellite (should we split this to LEO objects?)
        artefacts_all@radon.Streak; % everything else that is not astronomical
        too_loongs_all@radon.Streak; % LEO objects and very long artefacts
        
    end
    
    properties % inputs/outputs
        
        filenames; % cell of strings
        
        image;
        image_subtracted; % full frame, non filtered image, with streaks removed
        image_processed; % current image or section that undergoes processing and streak subtractions
        image_output; % image after truncating bright pixels and removing streaks
        image_sections = {}; % cell used to store sections of the image
        offset_sections = {};
        
        input_var;
        var_sections = {};
        
        runtime_preprocess;
        runtime_preprocess_total;
        runtime_streaks;
        runtime_streaks_total;
        
    end
    
    properties % switches/controls
        
        use_background_subtraction = 1;
        minimal_value = -300;
        use_trim = 1;
        trim_pixels = 30;
        final_size = 3072;
        
        threshold = 10;
        step_base = 2;
        
        use_scan_thresh = 1;
        use_conv = 1;
        iter_max = 5; % can increase this I think
        
        use_sections = 0;
        section_size = 1024;
        
        use_keep_too_longs = 1;
        use_keep_artefacts = 1;
        use_keep_satellites = 1;
        use_keep_asteroids = 1;
        use_keep_images = 0;
        
        debug_bit = 1;
        
        brake_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        input_psf;
        psf;
        psf_sigma;
        var_scalar;
        list;
        
    end
    
    properties(Hidden=true)
       
        current_filename = '';
        frame_index = 1;
        section_index = 1;
        
        total_frames;
        total_sections;
        
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
            obj.too_loongs_all = radon.Streak.empty;

            obj.clear;
            
            obj.frame_index = 1;
            
            obj.runtime_preprocess_total = [];
            obj.runtime_streaks_total = [];
            
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
            obj.too_longs = radon.Streak.empty;
            
            obj.runtime_preprocess = [];
            obj.runtime_streaks = [];
            
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
        
        function val = get.list(obj) % shorthand 
            
            val = obj.streaks_all;
            
        end
        
        function val = current_threshold(obj)
            
            val = obj.finder.threshold;
            
        end
        
    end
    
    methods % setters
        
        function set.input_psf(obj, val)
            
            obj.finder.input_psf = val;
            
        end
        
        function set.image(obj, val)
            
            obj.image = val;
%             obj.image_subtracted = val;
            
        end
        
    end
    
    methods % calculations
        
        function run(obj, varargin)
            
            input = util.text.InputVars;
            input.input_var('filenames', {}, 'files', 'names');
            input.input_var('reset', false);
            input.scan_vars(varargin{:});
            
            if ~isempty(input.filenames) && iscellstr(input.filenames)
                obj.filenames = input.filenames;
            end
            
            if input.reset
                obj.reset;
            end
            
            obj.brake_bit = 0;
            
            obj.total_frames = length(obj.filenames);
            
            for ii = obj.frame_index:length(obj.filenames)
                
                if obj.brake_bit
                    return;
                end
                
                if ~isempty(obj.gui), obj.gui.update; end
                
                [~, name, ext] = fileparts(obj.filenames{ii});
                name = [name, ext];
                
                t = tic;
                
                I = fitsread(obj.filenames{ii}, 'image');
                info = fitsinfo(obj.filenames{ii});
                seeing = info.Image.Keywords{strcmpi(info.Image.Keywords(:,1), 'seeing'),2};
                psf_sigma = seeing./2.355;
                
                I(I<obj.minimal_value) = NaN;
                
                if obj.use_trim
                    I = util.img.crop2size(I, obj.final_size-obj.trim_pixels);
                end
                
                I = util.img.pad2size(util.img.crop2size(I, obj.final_size), obj.final_size, NaN);
                
                [M,V] = util.img.tile_stats(I, 'tile', 32, 'varnan', 0);
                
                if obj.use_background_subtraction
                    I = I - M;
                end
                
                obj.runtime_preprocess = toc(t);
                obj.runtime_preprocess_total = obj.runtime_preprocess_total + toc(t);
                
                if obj.debug_bit, fprintf('Time to load and preprocess frame %d was %f seconds.\n', ii, obj.runtime_preprocess); end
                
                obj.input('image', I, 'psf', psf_sigma, 'variance', V, 'frame', ii, 'filename', name);
                
            end
            
            
        end
        
        function input(obj, varargin)
            
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('image', []);
            input.input_var('variance', [], 'input_variance');
            input.input_var('psf', [], 'input_psf');
            input.input_var('filename', '');
            input.input_var('batch_num', []);
            input.input_var('frame_num', []);
            input.input_var('show_bit', 0);
            input.scan_vars(varargin{:});
            
            if isempty(input.image)
                return; % quietly quit if there are no images
            end
            
            obj.current_filename = input.filename;
            obj.frame_index = input.frame_num;
            
            t = tic;
            
            obj.clear;
            
            obj.image = input.image;
            obj.image_subtracted = input.image;
            
            if ~isempty(input.psf)
                obj.input_psf = input.psf;
            end
            
            if ~isempty(input.variance)
                obj.input_var = input.variance;
            end
            
            if obj.use_sections
                
                obj.makeSections;
                
                obj.total_sections = length(obj.image_sections);
                
                for ii = 1:length(obj.image_sections)
                    
                    obj.image_processed = obj.preprocess(obj.image_sections{ii});
                    obj.runSection('offset', obj.offset_sections{ii}, 'section', ii,...
                        'variance', obj.var_sections{ii}, 'psf', input.psf, ...
                        'filename', input.filename, 'batch', input.batch_num, 'frame', input.frame_num, 'show', input.show_bit);
                    
                end
                
                if obj.gui.check
                    obj.showOriginal('ax', obj.gui.axes_image, 'bias', 0, 'dyn', obj.threshold);
                    drawnow;
                end
                
            else
                
                obj.image_processed = obj.preprocess(obj.image);
                obj.runSection('variance', obj.input_var, 'psf', input.psf,...
                    'filename', input.filename, 'batch', input.batch_num, 'frame', input.frame_num, 'show', input.show_bit);
                
            end
            
            for kk = 1:length(obj.streaks)
            
                if (obj.streaks(kk).is_too_long && obj.use_keep_too_longs) ||...
                            (obj.streaks(kk).is_artefact && obj.use_keep_artefacts) ||...
                            (obj.streaks(kk).is_satellite && obj.use_keep_satellites) ||...
                            (obj.streaks(kk).is_asteroid && obj.use_keep_asteroids) 
                        
                    obj.streaks_all = [obj.streaks_all obj.streaks(kk)];
            
                end
            
                if obj.streaks(kk).is_too_long && obj.use_keep_too_longs
                    obj.too_longs_all = [obj.too_longs_all obj.streaks(kk)];                
                end
                
                if obj.streaks(kk).is_artefact && obj.use_keep_artefacts
                    obj.artefacts_all = [obj.artefacts_all obj.streaks(kk)];                
                end
                
                if obj.streaks(kk).is_satellite && obj.use_keep_satellites
                    obj.satellites_all = [obj.satellites_all obj.streaks(kk)];
                end
                
                if obj.streaks(kk).is_asteroid && obj.use_keep_asteroids
                    obj.asteroids_all = [obj.asteroids_all obj.streaks(kk)];
                end
            
            end
            
            obj.runtime_streaks = toc(t);
            obj.runtime_streaks_total = obj.runtime_streaks_total + toc(t);
            
            if obj.debug_bit, fprintf('Time to find streaks in frame %d was %f seconds\n', obj.frame_index, obj.runtime_streaks); end
            
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
                        obj.var_sections{end+1} = obj.input_var(axis1(ii):end_point1, axis2(jj):end_point2);
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
            
            import util.stat.sum2;
            
            if ~isempty(obj.gui), obj.gui.update; end
            
            input = util.text.InputVars;
            input.input_var('variance',[]);
            input.input_var('filename', '');
            input.input_var('batch_num', []);
            input.input_var('frame_num', []);
            input.input_var('section', 1, 'section_number');
            input.input_var('offset', [0 0]);
            input.input_var('show_bit', 0);
            input.scan_vars(varargin{:});
            
            obj.section_index = input.section;
            
            I = obj.image_processed; % filtered image! 
            V = input.variance; % scalar value of the variance (or median of the var map)
            
            if isempty(V)
                V = 1; 
            elseif ismatrix(V)
                V = util.stat.median2(V);
            end
                        
            if obj.use_scan_thresh
                log_values = log2(obj.threshold)./log2(obj.step_base):log2(util.stat.max2(I)./sqrt(V.*sum2(obj.psf)))./log2(obj.step_base);
                if isempty(log_values)
                    thresh_values = obj.threshold;
                else
                    thresh_values = flip(obj.step_base.^log_values);
                    thresh_values = [thresh_values(1)*obj.step_base thresh_values]; % add one more level for bright streaks in uniform images
                end
            else
                thresh_values = obj.threshold;
            end
            
            for ii = 1:length(thresh_values)
                
                if obj.gui.check
                    obj.showProcessed('thresh', thresh_values(ii), 'section', input.section);
                    drawnow;
                end
                
                obj.finder.threshold = thresh_values(ii);
                
                mask = (I./sqrt(V.*sum2(obj.psf)))>thresh_values(ii);
                
                I(mask) = thresh_values(ii).*sqrt(V.*sum2(obj.psf));
                
                obj.image_output = I;
                
                if obj.debug_bit>1, fprintf('  Threshold= %4.2f | mask= % 6d (%8.6f%%)\n', thresh_values(ii), nnz(mask), nnz(mask)./numel(I)*100); end 
                
                temp_streaks = radon.Streak.empty;

                if ~isempty(obj.gui)
                    obj.gui.update;
                    drawnow;
                end
                
                for jj = 1:obj.iter_max
                    
                    new_streaks = obj.finder.input('image', I, 'variance', input.variance,...
                         'was_convolved', obj.use_conv, 'original_image', obj.image, ...
                         'filename', input.filename, 'batch_num', input.batch_num,...
                         'frame_num', input.frame_num, 'section_num', input.section, 'offset', input.offset);
                    
                    if isempty(new_streaks)
                        break;
                    end
                    
                    for kk = 1:length(new_streaks)
                        
                        new_streaks(kk).prof.input(obj.image_subtracted, 'show', input.show_bit);
                        obj.class.input(new_streaks(kk));
                        
                    end
                    
                    assert(length(new_streaks)<=2, 'The following code is only suitable for getting <=2 streaks from finder');
                    
                    if length(new_streaks)==2 
                        
                        s1 = new_streaks(1);
                        s2 = new_streaks(2);
                        
                        if abs(s1.th-s2.th)<5 || ...
                                sqrt((s1.midpoint_x-s2.midpoint_x).^2+(s1.midpoint_y-s2.midpoint_y).^2)<40
                            new_streaks = s1; % is there a better way to choose one of these streaks? maybe the updated photometry S/N
                        
                        end
                        
                    end
                        
                    for kk = 1:length(new_streaks)
                        
                        I = new_streaks(kk).subtractStreak('image', I, 'replace', NaN); % used as image_output (filtered and subtracted)
                        obj.image_subtracted = new_streaks(kk).subtractStreak(... % full image, unfiltered, subtracted
                            'image', obj.image_subtracted, 'replace', NaN, 'offset', new_streaks(kk).offset_original);
                        
                        if obj.use_keep_images==0
                            new_streaks(kk).image = [];
                            new_streaks(kk).original_image = [];
                            new_streaks(kk).radon_image = [];
                            new_streaks(kk).subframe = [];
                        end
                        
                        % add the new streaks to the list
                        
                        obj.streaks = [obj.streaks new_streaks(kk)];
                        
                        if new_streaks(kk).is_too_long && obj.use_keep_too_longs
                            obj.too_longs = [obj.too_longs new_streaks(kk)];
                        end
                        
                        if new_streaks(kk).is_artefact && obj.use_keep_artefacts
                            obj.artefacts = [obj.artefacts new_streaks(kk)];
                        end
                        
                        if new_streaks(kk).is_satellite && obj.use_keep_satellites
                            obj.satellites = [obj.satellites new_streaks(kk)];
                        end
                        
                        if new_streaks(kk).is_asteroid && obj.use_keep_asteroids
                            obj.asteroids = [obj.asteroids new_streaks(kk)];
                        end
                        
                        if ~isempty(obj.second_fig) && isvalid(obj.second_fig) && isa(obj.second_fig, 'matlab.ui.Figure')
                            if ~new_streaks(kk).is_too_long
                                new_streaks(kk).showAll('parent', obj.second_fig);
                                drawnow;
                            end
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
                        fprintf('Found streaks in frame %d: SNR= %s | thresh= %s | L= %s\n', obj.finder.frame_num, ...
                            util.text.print_vec([temp_streaks.snr]), util.text.print_vec([temp_streaks.threshold]),...
                            util.text.print_vec([temp_streaks.L])); 
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
%                 obj.streaks(ii).drawGuidelines('ax', input.axes);
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

