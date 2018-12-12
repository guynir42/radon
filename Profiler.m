classdef Profiler < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        streak@radon.Streak; % point back to the streak that contains this profiler
        
    end
    
    properties % inputs/outputs
        
        image;
        expT = 30; % default value
        
        image_filt;
        image_rot;
        image_crop;
        
        trunc_mask; % area of the filtered image that need to be truncated
        star_mask; % area of the filtered image that has stars in it
        
        trunc_mask_rot;
        star_mask_rot;
        trunc_mask_crop;
        star_mask_crop;
        
%         image_bin;
        
        width;
        bin_size;
        
        filter_width;
        filter_snr;
        filter_loss;
        filter_edge;
        
        angle_correction;
        
        amplitudes;
        amp_lowpass;
        offsets;
        spreads;
        line_start;
        line_end;
        variability;

        glint_power;
        
        lateral_profile;
        lateral_spread;
        
        fit_period;
        fit_amplitude;
        fit_slope;
        fit_residual;
                
        type = '';
        
    end
    
    properties % switches/controls
        
        max_length = 500;
        min_snr = 50;
        nsigma_width = 10;
        
        lowpass_freq = 0.1;
        
        debug_bit = 1;
        show_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
        
        lateral_profile_regular;
        lateral_profile_filter;
        lateral_profile_overfilter;
        
        fit_options;
        fit_type;
        
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Profiler(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'radon.Profiler')
                if obj.debug_bit, fprintf('Profiler copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            elseif ~isempty(varargin) && isa(varargin{1}, 'radon.Streak')
%                 if obj.debug_bit, fprintf('Profiler streak constructor v%4.2f\n', obj.version); end
                obj.streak = varargin{1};
                obj.initialize;
            else
                if obj.debug_bit, fprintf('Profiler constructor v%4.2f\n', obj.version); end
                obj.initialize;
            end
            
        end
        
    end
    
    methods % reset/clear
        
        function initialize(obj)
            
            obj.reset;
            
        end
        
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
            
            obj.filter_width = [];
            obj.filter_snr = [];
            obj.filter_loss = [];
            obj.filter_edge = [];
            
            obj.angle_correction = [];
            
            obj.amplitudes = [];
            obj.offsets = [];
            obj.spreads =[];
            obj.line_start = [];
            obj.line_end = [];
            obj.variability = [];

            obj.fit_period = [];
            obj.fit_amplitude = [];
            obj.fit_slope = [];
            obj.fit_residual = [];
            
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
        
        function input(obj, varargin)
           
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('image', []);
            input.scan_vars(varargin);
                        
            if isempty(obj.streak)
                error('Profiler must have a Streak object to run!');
            end
                        
            if isempty(input.image)
                input.image = obj.streak.original_image;
            end
            
            if isempty(input.image) 
                error('Streak does not have an original image and no input image was given!');
            end
            
            obj.clear;
            
            obj.image = obj.streak.getCutout('width', obj.streak.radon_dx.*2, 'image', input.image);
%             obj.image = util.img.pad2size(obj.image, 2*ceil(size(obj.image)/2)); % round image dimensions to even sizes
            
            obj.width = round(obj.streak.psf_sigma*2*obj.nsigma_width+1);
            obj.bin_size = round(obj.streak.psf_sigma*2); % do we even need this anymore?            
            if obj.bin_size<1, obj.bin_size = 1; end
            
            % test number 1: length of streak from FRT (don't continue
            % testing very long streaks aka LEO satellites)
            if obj.streak.L>obj.max_length
                % need to do anything else?
                return;
            end
            
            obj.calcFiltering; % test 2: best PSF sigma to filter with
            obj.calcBestAngle; % test 3: find the best rotation angle around the FRT output
            obj.calcLength; % test 4: find the real length of the streak
            obj.calcAmpOffsetSpread; % test 5: find the amplitudes, offsets and spreads of the streak
            obj.calcVariability;
            obj.calcLateralSpread;
            obj.calcPeriodicity; % test 6: find the best periodicity in the signal
            
%                 obj.fitLine;
            
        end
        
        function calcFiltering(obj)
            
            import util.stat.sum2;
            import util.stat.max2;
            
            filter_range = -3.5:0.25:3.5;
            sig = filter_range + obj.streak.psf_sigma;
            sig = sig(sig>=0.3);
            
            RV = radon.frt(ones(size(obj.image)).*obj.streak.noise_var); % radon space variance (for normalization)
            
            I = obj.image;
            I(obj.star_mask) = NaN;
            I = I - util.stat.median2(I);
            I(isnan(I)) = 0;
            
            obj.image = I;
            
            SNR = zeros(1,length(sig));
            
            for ii = 1:length(sig)
                
                P = util.img.gaussian2(sig(ii), 'norm', 2);
                
                If = filter2(P,I);
                If(If./sqrt(obj.streak.noise_var)>obj.streak.threshold) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
                
                if obj.show_bit
                    util.plot.show(If, 'bias', 0);
                    pause(0.5);
                end
                
                R = radon.frt(If, 'transpose', 0); 
                RT = radon.frt(If, 'transpose', 1); 
                
                SNR(ii) = max(max2(RT./sqrt(RV)), max2(R./sqrt(RV)))./sqrt(sum2(P));
                
                if sig(ii)==obj.streak.psf_sigma
                    expected_snr = SNR(ii);
                end
                
            end
            
            obj.filter_width = sig;
            obj.filter_snr = SNR;
            
            [mx, idx] = max(SNR);
            obj.filter_loss = abs(mx-expected_snr)./mx; 
            if idx==1 || idx==length(SNR)
                obj.filter_edge = 1;
            else
                obj.filter_edge = 0;
            end
            
            if obj.show_bit
                plot(sig, SNR);
                xlabel('PSF sigma');
                ylabel('relative SNR');
                pause(2);
            end
            
        end
        
        function calcBestAngle(obj)
            
            % can later make these inputs to the profiler
            range = 3;
            step = 0.5;
            
            angles = -range:step:range;
            
            score = zeros(length(angles),1);
            pos = zeros(length(angles),1);
            
            obj.image_filt = filter2(obj.streak.psf, obj.image);
            obj.trunc_mask = obj.image_filt./sqrt(obj.streak.noise_var)>obj.streak.threshold;
            obj.star_mask = imdilate(obj.image_filt./sqrt(obj.streak.noise_var)>obj.streak.threshold*2, ones(9));
            
%             obj.image_filt(obj.trunc_mask) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
            I = obj.image_filt;
            I(obj.trunc_mask) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
            I(obj.star_mask) = 0;

            for ii = 1:length(angles)
                
                Ir = imrotate(I, angles(ii)+obj.streak.th, 'bilinear','crop');
                [score(ii), pos(ii)] = max(sum(Ir,2));
                
                if obj.show_bit                    
                    plot(sum(Ir,2));
                    hold on;
                    xlabel(['angle= ' num2str(angles(ii)) ' | score= ' num2str(score(ii))]);
                    pause(0.5);
                end
                
            end
            
            if obj.show_bit
                hold off; 
                plot(angles, score);
                pause(2);
            end
            
            [~,idx] = max(score);
            pos = round(median(pos)); % where the streak on the y axis
            
            obj.angle_correction = angles(idx);
            
            % make a rotated, then a cropped image
            obj.image_rot = imrotate(obj.image, obj.streak.th+obj.angle_correction, 'bilinear','crop');
            obj.trunc_mask_rot = filter2(obj.streak.psf, obj.image_rot)./sqrt(obj.streak.noise_var)>obj.streak.threshold;
            obj.star_mask_rot = imdilate(filter2(obj.streak.psf, obj.image_rot)./sqrt(obj.streak.noise_var)>obj.streak.threshold*2, ones(9));
            
            top = pos - floor(obj.width/2-1);
            if top<1, top = 1; end
            
            bottom = pos + floor(obj.width/2+1);
            if bottom>size(obj.image_rot,1), bottom = size(obj.image_rot,1); end
            
            obj.image_crop = obj.image_rot(top:bottom,:);
            obj.trunc_mask_crop = obj.trunc_mask_rot(top:bottom,:);
            obj.star_mask_crop = obj.star_mask_rot(top:bottom,:);
            
        end
        
        function calcLength(obj)
            
            % if we want to used some kind of binning before length estimates
%             obj.image_bin = filter2(ones(1,obj.bin_size), obj.image_crop);
%             obj.image_bin = obj.image_bin(:,1:obj.bin_size:end);
%             I = filter2(sum(obj.streak.psf,2), obj.image_bin); % binned along the streak, matched-filtered in the lateral direction...
%             obj.image_filt = filter2(obj.streak.psf, obj.image_rot);

            I = obj.image_crop;
            If = filter2(obj.streak.psf, I);
            If(obj.trunc_mask_crop) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
            If(obj.star_mask_crop) = 0;
            
            N = size(I,2);
            
            best_start = 1;
            best_end = N;
            best_snr = 0;
            
            for ii = 1:N % we can improve this a little bit using smart additions (is this slowing us down?)
                for jj = ii:N
                    
                    noise = sqrt(jj-ii+1); % how many pixels in the line
                    signal = sum(max(If(:, ii:jj)), 'omitnan');
                    snr = signal./noise;
                    
                    if snr>best_snr
                        best_snr = snr;
                        best_start = ii;
                        best_end = jj;
                    end
                    
                end
            end
            
%             obj.line_start = best_start.*obj.bin_size;
%             obj.line_end = best_end.*obj.bin_size;
            
            obj.line_start = best_start;
            obj.line_end = best_end;
            
        end
        
        function calcAmpOffsetSpread(obj)
            
            I = obj.image_crop;
            I(obj.star_mask_crop) = 0;
            
            for ii = 1:size(obj.image_crop,2)

                slice = I(:,ii);
                y = (1:length(slice))'-floor(length(slice)/2)-1;
                obj.amplitudes(ii) = max(slice);
                
                slice(slice<0) = 0;
                S = sum(slice);
                if S==0, S = 1; end % this shouldn't happen but if it does, I prefer to get offsets and spreads equal zero than Inf
                
                obj.offsets(ii) = sum(slice.*y, 'omitnan')./S;
                obj.spreads(ii) = sqrt(sum(slice.*(y-obj.offsets(ii)).^2, 'omitnan')./S);
                
            end
                        
            pixel_width = mean(obj.spreads(obj.line_start:obj.line_end));
            new_value = 2*pixel_width./obj.streak.psf_sigma;
            
            if new_value>obj.streak.subtract_psf_widths
                obj.streak.subtract_psf_widths = new_value;
            end
            
        end
        
        function calcVariability(obj)
            
            a = obj.amplitudes(obj.line_start:obj.line_end);
            
            v = a + obj.streak.noise_var; % background and source noise together (assuming GAIN=1)
            
            b = a./sqrt(v); % normalize to unity noise
            
            obj.variability = std(b);
            
%             obj.amp_lowpass = lowpass(obj.amplitudes,obj.lowpass_freq);
            
            obj.amp_lowpass = lowpass(a, obj.lowpass_freq, 1);
%             obj.amp_lowpass = lowpass(obj.amplitudes, obj.lowpass_freq, 1);
%             obj.amp_lowpass = obj.amp_lowpass(obj.line_start:obj.line_end);

            obj.glint_power = std(obj.amp_lowpass)./sqrt(median(abs(obj.amp_lowpass)));
            
        end
        
        function calcLateralSpread(obj)
            
%             I = filter2(obj.streak.psf, obj.image_crop);
%             I(obj.trunc_mask_crop) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
            
            I = obj.image_crop(:,obj.line_start:obj.line_end);
            slice = sum(I,2);

            slice(slice<0) = 0;
            
            y = (1:length(slice))'-floor(length(slice)/2)-1;
            
            m1 = sum(slice.*y)./sum(slice);
            m2 = sum(slice.*(y-m1).^2)./sum(slice);
            
            obj.lateral_spread = sqrt(m2);
            obj.lateral_profile = slice;
            
        end
        
        function calcPeriodicity(obj)
            
            values = obj.amplitudes(obj.line_start:obj.line_end)'./sqrt(obj.streak.noise_var);
            times = linspace(0,obj.expT,length(values))';
            dt = times(2)-times(1);
            
            func = @(p) util.stat.lscov_sinus(times, values, 1./p);
            
%             p_found = fminsearch(func,3);

            periods = 5*dt:dt:obj.expT*10;
            res = zeros(length(periods),1);
            
            for ii = 1:length(periods)
                res(ii) = func(periods(ii));
            end
            
            [min_res, idx] = min(res);
            p_found = periods(idx);
            
            if obj.show_bit
                plot(periods, res);
                xlabel('period [seconds]');
                ylabel('residuals');            
                title(['best period= ' num2str(p_found) ' | res= ' num2str(min_res)]);
                pause(2);
            end
            
            [residual, pars, model_data] = func(p_found);
            
            if obj.show_bit
                plot(times,values, '.', times, model_data);
                legend({'data', 'fit'});
                xlabel('time [sec]');
                ylabel('brightness [normalized]');
                title(['period= ' num2str(p_found) ' | res= ' num2str(residual)]);
                pause(2);
            end
            
            obj.fit_period = p_found;
            obj.fit_amplitude = sqrt((pars(3).^2+pars(4).^2)*obj.streak.noise_var);
            obj.fit_slope = pars(2);
            obj.fit_residual = residual;
            
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
        
        function calcFiltering2(obj) % to be depricated! 
            
            import util.stat.sum2;
            
            % filter the image using 
            if isempty(obj.streak.original_image) && obj.streak.was_convolved
                P = obj.streak.psf;
                obj.image_filt = obj.image_rot; % this means the streak only has the processed, filtered, image
                P_over = util.img.gaussian2(3.*obj.streak.psf_sigma, 'norm', 2);
                obj.image_overfilt = filter2(P_over, obj.image_rot); % use slightly sharper PSF for overfiltering
            else
                P = obj.streak.psf;
                obj.image_filt = filter2(P, obj.image_rot);
                P_over = util.img.gaussian2(4.*obj.streak.psf_sigma, 'norm', 2);
                obj.image_overfilt = filter2(P_over, obj.image_rot);
            end
            
%             correction = sqrt(obj.streak.L.*obj.streak.noise_var.*max(abs(cosd(obj.streak.th)), abs(sind(obj.streak.th))));
            correction = sqrt(obj.streak.L.*obj.streak.noise_var);
            
            obj.lateral_profile_filter = sum(obj.image_filt,2)./sqrt(sum2(P))./correction;
            obj.streak.filtered_snr = max(obj.lateral_profile_filter);
            
            obj.lateral_profile_overfilter = sum(obj.image_overfilt,2)./sqrt(sum2(P_over))./correction;
            obj.streak.overfiltered_snr = max(obj.lateral_profile_overfilter);
            
            if isempty(obj.streak.original_image) && obj.streak.was_convolved
                obj.lateral_profile_regular = obj.lateral_profile_filter;
                obj.streak.unfiltered_snr = obj.streak.filtered_snr;
            else
                obj.lateral_profile_regular = sum(obj.image_rot,2)./correction;
                obj.streak.unfiltered_snr = max(obj.lateral_profile_regular);
            end
            
            if obj.show_bit
            
                plot(obj.lateral_profile_regular);
                hold on;
                plot(obj.lateral_profile_filter);
                plot(obj.lateral_profile_overfilter);
                hold off;                
                legend({'Non-filtered', 'Filtered', 'Over-filtered'});
%                 xlabel(sprintf('REG: mx= %4.2f | std= %4.2f --- FILT: mx= %4.2f | std= %4.2f', ...
%                     score_regular.*std_regular, std_regular, score_filter.*std_filter, std_filter));
            
                pause(2);
                
            end
            
        end
        
        function [score, sd] = calcLateralProfileScore(obj, prof) % to be depricated! this gives a score that is much too low for the filtered image (noise is overestimated...)
            
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
        
    end
    
    methods % plotting tools / GUI
        
        function plotFiltering(obj, varargin)
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            plot(input.axes, obj.filter_width, obj.filter_snr);
            xlabel(input.axes, 'Filter \sigma');
            ylabel(input.axes, 'S/N'); 
            
            util.plot.inner_title(['loss= ' num2str(obj.filter_loss)], 'Position', 'Bottom', 'axes', input.axes);
            
        end
        
        function plotAmplitudes(obj, varargin)
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            prev_hold_state = input.axes.NextPlot;
            plot(input.axes, 1:length(obj.amplitudes), obj.amplitudes);
            hold(input.axes, 'on');
            plot(input.axes, obj.line_start:obj.line_end, obj.amplitudes(obj.line_start:obj.line_end));
%             plot(input.axes, obj.line_start:obj.line_end, lowpass(obj.amplitudes(obj.line_start:obj.line_end), obj.lowpass_freq));
            plot(input.axes, obj.line_start:obj.line_end, obj.amp_lowpass);
            input.axes.NextPlot = prev_hold_state;
            
            xlabel(input.axes, 'position along streak');
            ylabel(input.axes, 'amplitude');
            
            util.plot.inner_title(['var= ' num2str(obj.variability)], 'Position', 'Bottom', 'axes', input.axes);
            util.plot.inner_title(['glint= ' num2str(obj.glint_power)], 'Position', 'Top', 'axes', input.axes);
            
        end
        
        function plotLateralProfile(obj, varargin)
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            plot(input.axes, obj.lateral_profile);
            
            xlabel(input.axes, 'lateral position');
            ylabel(input.axes, 'summed intensity');
            
            util.plot.inner_title(['spread= ' num2str(obj.lateral_spread)], 'Position', 'Bottom', 'axes', input.axes);
            
        end
        
        function plotSpreads(obj, varargin)
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
%             prev_hold_state = input.axes.NextPlot;
            x = obj.line_start:obj.line_end;
            y = obj.spreads(obj.line_start:obj.line_end);
            plot(input.axes, x, y);
%             hold(input.axes, 'on');
            
%             plot(input.axes, obj.line_start:obj.line_end, obj.offsets(obj.line_start:obj.line_end));
%             input.axes.NextPlot = prev_hold_state;
            
            input.axes.YLim = [0 1.2*max(y)];

            xlabel(input.axes, 'position within streak');
            ylabel(input.axes, 'spread/offset');
            
            util.plot.inner_title(['mean spread= ' num2str(mean(y))], 'Position', 'Bottom', 'axes', input.axes);
                        
        end
        
    end    
    
end

