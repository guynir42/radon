classdef Profiler < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        streak@radon.Streak;
        
    end
    
    properties % inputs/outputs
        
        image;
        image_rot;
        image_crop;
        image_filt;
        image_bin;
        width;
        bin_size;
        
        line_start;
        line_end;
        
        amplitudes;
        offsets;
        spreads;
        
        type = '';
        
    end
    
    properties % switches/controls
        
        max_length = 500;
        min_snr = 50;
        nsigma_width = 10;
        
        debug_bit = 1;
        show_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        lateral_profile_regular;
        lateral_profile_filter;
        
        fit_options;
        fit_type;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Profiler(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'radon.Profiler')
                if obj.debug_bit, fprintf('Profiler copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Profiler constructor v%4.2f\n', obj.version); end
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function reset(obj)
            
            obj.clear;
            
        end
        
        function clear(obj)
            
            obj.image = [];
            obj.image_rot = [];
            obj.image_crop = [];
            obj.image_filt = [];
            obj.width = [];
            obj.bin_size = [];
            
            obj.line_start = [];
            obj.line_end = [];
            
            obj.amplitudes = [];
            obj.offsets = [];
            obj.spreads =[];
            
            obj.fit_options = [];
            obj.fit_type = [];
            
            obj.lateral_profile_regular = [];
            obj.lateral_profile_filter = [];
            
        end
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, streak)
           
            if nargin==0 || isempty(streak) || ~isa(streak, 'radon.Streak')
                error('Must give a "radon.Streak" object to "input"');
            end
            
            obj.clear;
            
            obj.streak = streak;
            obj.image = obj.streak.getCutout;
            obj.width = round(obj.streak.psf_sigma*2*obj.nsigma_width+1);
            obj.bin_size = round(obj.streak.psf_sigma*2);
            if obj.bin_size<1, obj.bin_size = 1; end
            
            obj.image_rot = util.fft.RotateImage(util.img.pad2size(obj.image, 2*ceil(size(obj.image)/2)), obj.streak.th);
            
            obj.image_crop = util.img.crop2size(obj.image_rot, [obj.width, size(obj.image_rot,2)]);
            
            obj.image_filt = filter2(obj.streak.psf, obj.image_rot);
            
            if obj.testFiltering
                
                if obj.debug_bit, disp('Calculating fit to streak segments'); end
                
                obj.image_bin = filter2(ones(1,obj.bin_size), obj.image_crop);
                obj.image_bin = obj.image_bin(:,1:obj.bin_size:end);
                obj.findLength;
                obj.estimateIntensityWidth;
%                 obj.fitLine;
            
                obj.streak.amplitdues = obj.amplitudes;
                obj.streak.offsets = obj.offsets;
                obj.streak.spreads = obj.spreads;
                obj.streak.line_start = obj.line_start;
                obj.streak.line_end = obj.line_end;
                pixel_width = mean(obj.spreads(obj.line_start:obj.line_end));
                obj.streak.subtract_psf_widths = 2*pixel_width./obj.streak.psf_sigma;
                obj.streak.is_astronomical = 1;
                
%                 if pixel_width>2.*obj.streak.psf_sigma
%                     obj.streak.is_astronomical = 0;
%                 end
                
            else
                if obj.debug_bit, disp('Streak marked as not-astronomical. Skipping line fitting'); end
                obj.type = 'non-astronomical';
                obj.streak.is_astronomical = 0;
            end
            
        end
        
        function is_astronomical = testFiltering(obj)
            
            obj.lateral_profile_regular = sum(obj.image_rot,2);
            score_regular = max(obj.lateral_profile_regular);
%             [score_regular, std_regular] = obj.calcLateralProfileScore(obj.lateral_profile_regular);
            
            obj.lateral_profile_filter = sum(obj.image_filt,2);
            score_filter = max(obj.lateral_profile_filter);
%             [score_filter, std_filter] = obj.calcLateralProfileScore(obj.lateral_profile_filter);
            
            if score_filter - sqrt(score_filter) > score_regular
                is_astronomical = 1;
            else
                is_astronomical = 0;
            end
            
            if obj.show_bit
            
                plot(obj.lateral_profile_regular);
                hold on;
                plot(obj.lateral_profile_filter); 
                hold off;                
                legend({'Non-filtered', 'Filtered'});
%                 xlabel(sprintf('REG: mx= %4.2f | std= %4.2f --- FILT: mx= %4.2f | std= %4.2f', ...
%                     score_regular.*std_regular, std_regular, score_filter.*std_filter, std_filter));
            
                pause(2);
                
            end
            
        end
        
        function [score, sd] = calcLateralProfileScore(obj, prof) % this gives a score that is much too low for the filtered image (noise is overestimated...)
            
            prof = util.vec.torow(prof);
            
            [mx, idx] = max(prof);
            prof_noise = [];
            if idx>obj.width
                prof_noise = [prof_noise prof(1:idx-obj.width)];
            end
            
            if idx<=length(prof)-obj.width
                prof_noise = [prof_noise prof(idx+obj.width:end)];
            end
            
            sd = std(prof_noise);
            
            score = mx./sd;
            
        end
        
        function findLength(obj)
            disp('findLength');
            I = filter2(sum(obj.streak.psf,2), obj.image_bin); % binned along the streak, matched-filtered in the lateral direction...
            N = size(I,2);
            center = floor(size(I,1)./2)+1;
            
            best_start = 1;
            best_end = N;
            best_snr = 0;
            
            for ii = 1:N
                for jj = ii:N
                    
                    noise = sqrt(jj-ii+1); % how many pixels in the line
                    signal = sum(I(center, ii:jj));
                    snr = signal./noise;
                    
                    if snr>best_snr
                        best_snr = snr;
                        best_start = ii;
                        best_end = jj;
                    end
                    
                end
            end
            
            obj.line_start = best_start.*obj.bin_size;
            obj.line_end = best_end.*obj.bin_size;
            
        end
        
        function estimateIntensityWidth(obj)
            disp('estimateIntensityWidth');
%             for ii = obj.line_start:obj.line_end
            for ii = 1:size(obj.image_crop,2)

                slice = (obj.image_crop(:,ii));
                y = (1:length(slice))'-floor(length(slice)/2);
                obj.amplitudes(ii) = max(slice);
                obj.offsets(ii) = sum(abs(slice).*y, 'omitnan')./sum(abs(slice), 'omitnan');
                obj.spreads(ii) = sqrt(sum(abs(slice).*y.^2, 'omitnan')./sum(abs(slice), 'omitnan'));
                
            end
            
        end
        
        function fitLine(obj)
            
            obj.fit_type = fittype('c + a*exp(-0.5*(x-b)^2/sig^2)', 'Problem', 'sig', 'Independent', 'x', 'Coefficients', {'a', 'b', 'c'});
            
            obj.fit_options = fitoptions('Method', 'NonLinearLeastSquares', 'Robust', 'Bisquare', 'Display', 'off', ...
                'Lower', [0,-obj.width*1/4,-Inf], 'Upper', [Inf, obj.width*1/4, Inf], 'StartPoint', [1,0,0]);
            
            I = obj.image_bin; 
            N = size(I,2);
            
            obj.amplitudes = zeros(1,N);
            obj.offsets = zeros(1,N);
            
            x = (1:obj.width)'-obj.width/2;
            
            for ii = 1:N
                
                fr = fit(x, I(:,ii), obj.fit_type, obj.fit_options, 'problem', obj.streak.psf_sigma);
                obj.amplitudes(ii) = fr.a;
                obj.offsets(ii) = fr.b;
                
                if obj.show_bit
                   
                    plot(fr, x, I(:,ii));
                    util.plot.inner_title(['ii= ' num2str(ii)], 'Position', 'Bottom');
                    drawnow;
                    
                end
                
            end
            
        end
        
    end
    
    methods % plotting tools / GUI
        
    end    
    
end

