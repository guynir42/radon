classdef Profiler < handle

    properties(Transient=true)
        
    end
    
    properties % objects
        
        streak@radon.Streak; % point back to the streak that contains this profiler
        fitter@radon.Fitter;
        
    end
    
    properties % inputs/outputs
        
        image;
        im_offset; % what we need to translate coordinates in the cutout to the real 
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
        
        good_fitness;
        
        filter_width;
        filter_snr;
        filter_loss;
        filter_edge;
        filter_delta;
        
        edgyness;
        edge_map;
        edge_norm;
        
        angle_correction;
        
        amplitudes;
        amp_lowpass;
        offsets;
        spreads;
        line_start;
        line_end;
        variability;

        brightest_value;
        brightest_number;
        brightest_ratio;
        
        lateral_profile;
        lateral_skewness;
        lateral_spread;
        
        glint_power;
                
        runtime;
        
        % these are all going away soon...
        fit_period;
        fit_amplitude;
        fit_slope;
        fit_residual;
                
        type = '';
        
    end
    
    properties % switches/controls
        
        use_angle_cut = 1;
        
        use_star_mask = 1;
        use_trunc_mask = 1;
        
        max_length = 500;
        max_filter_delta = 0.5;
        max_filter_loss = 0.1;
        min_snr = 50;
        nsigma_width = 10;
        
        angle_step = 0.5;
        angle_range = 6.0;
        
        lowpass_freq = 0.08;
        
        debug_bit = 1;
        show_bit = 0;
        
    end
    
    properties(Dependent=true)
        
        l1; % shortcut property for "line_start" (also, defaults to 1)
        l2; % shortcut property for "line_end" (also, defaults to size of cutout)
        
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
            
            obj.good_fitness = [];
            
            obj.filter_width = [];
            obj.filter_snr = [];
            obj.filter_loss = [];
            obj.filter_edge = [];
            
            obj.edgyness = [];
            obj.angle_correction = [];
            
            obj.amplitudes = [];
            obj.offsets = [];
            obj.spreads =[];
            obj.line_start = [];
            obj.line_end = [];
            obj.glint_power = [];
            obj.variability = []; % maybe depricated soon

            obj.brightest_value = [];
            obj.brightest_number = [];
            obj.brightest_ratio = [];
            
            obj.lateral_profile = [];
            obj.lateral_skewness = [];
            obj.lateral_spread = [];
            
            obj.fit_options = [];
            obj.fit_type = [];
            
            % not sure if we need these...
            obj.fit_period = [];
            obj.fit_amplitude = [];
            obj.fit_slope = [];
            obj.fit_residual = [];
            
            obj.lateral_profile_regular = [];
            obj.lateral_profile_filter = [];
            
        end
        
    end
    
    methods % getters
        
        function val = get.l1(obj)
            
            if isempty(obj.line_start)
                val = 1;
            else 
                val = obj.line_start;
            end
            
        end
        
        function val = get.l2(obj)
            
            if isempty(obj.line_end)
                val = size(obj.image);
            else 
                val = obj.line_end;
            end
            
        end
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, varargin)
           
            input = util.text.InputVars;
            input.use_ordered_numeric = 1;
            input.input_var('image', []);
            input.input_var('offset', []);
            input.input_var('show_bit', 0);
            input.scan_vars(varargin);
            
            if isempty(obj.streak)
                error('Profiler must have a Streak object to run!');
            end
            
            t= tic;
            
            if isempty(input.image)
                if ~isempty(obj.streak.original_image)
                    input.image = obj.streak.original_image;
                else
                    input.image = obj.streak.image;
                    input.offset = [0 0];
                end
            end
            
            if isempty(input.offset)
                if all(size(input.image)==size(obj.streak.original_image))
                    input.offset = obj.streak.offset_original;
                else
                    input.offset = [0 0];
                end
            end
            
            if isempty(input.image) 
                error('Streak does not have an original image. Please input the original image...');
            end
            
            obj.clear;
            cleanup = onCleanup(@() obj.finishup(t));
            
            [obj.image, obj.im_offset] = obj.streak.getCutout('width', obj.streak.radon_dx.*4, 'image', input.image, 'offset', input.offset);
%             obj.image = util.img.pad2size(obj.image, 2*ceil(size(obj.image)/2)); % round image dimensions to even sizes
            obj.im_offset = obj.im_offset - obj.streak.offset_original;
            
            obj.width = round(obj.streak.psf_sigma*2*obj.nsigma_width+1);
            obj.bin_size = round(obj.streak.psf_sigma*2); % do we even need this anymore?            
            if obj.bin_size<1, obj.bin_size = 1; end
            
            % test number 1: length of streak from FRT (don't continue
            % testing very long streaks aka LEO satellites)
            if obj.streak.L>obj.max_length
                % need to do anything else?
                return;
            end
            
%             obj.calcFiltering; % best PSF sigma to filter with
%             obj.calcBestAngle; % find the best rotation angle around the FRT output
%             obj.calcLength; % find the real length of the streak
%             obj.calcBestAngle; % find the best rotation angle around the FRT output
%             obj.calcLength; % find the real length of the streak

            if obj.use_angle_cut && abs(mod(obj.streak.th,90))<1
                return;
            end

            obj.rotateCrop;
            obj.calcFiltering;
            
            if obj.filter_delta>obj.max_filter_delta || obj.filter_edge || obj.filter_loss>obj.max_filter_loss
                return;
            end
            
            obj.calcFit(input.show_bit);
            obj.rotateCrop;
            obj.calcAmpOffsetSpread; % find the amplitudes, offsets and spreads of the streak
%             obj.calcFiltering;
            obj.calcVariability;
            obj.calcLateralSpread;
            obj.calcBrightest;
            obj.calcEdgyness;
%             obj.calcPeriodicity; % test 6: find the best periodicity in the signal
            
%                 obj.fitLine;

%             obj.correctEdges;
            
        end
        
        function finishup(obj, t)
            
            obj.runtime = toc(t);
            
        end
        
        function calcFiltering(obj)
            
            import util.stat.sum2;
            import util.stat.max2;
            
            filter_range = -3.5:0.1:3.5;
            sig = filter_range + obj.streak.psf_sigma;
            sig = sig(sig>=0.3); 
            
            I = obj.image_crop;
            
            SNR = zeros(1,length(sig));
            
            for ii = 1:length(sig)
                
                P = util.img.gaussian2(sig(ii), 'norm', 2);
                
                If = filter2(P,I);
%                 If = util.img.truncate(If, obj.streak.threshold, obj.streak.noise_var);
                
%                 line = sum(If(:,obj.line_start:obj.line_end),2)./sqrt(obj.streak.noise_var.*(obj.line_end-obj.line_start).*sum2(P));
%                 [SNR(ii), ~] = max(line);
                SNR(ii) = sum(max(If(:,obj.line_start:obj.line_end)),2, 'omitnan')./sqrt(obj.streak.noise_var.*(obj.line_end-obj.line_start).*sum2(P));
                
            end
            
            [~,idx1] = min(abs(sig-obj.streak.psf_sigma_effective)); % the closest value to the "effective PSF" sigma
            expected_snr = SNR(idx1);
            
            obj.filter_width = sig;
            obj.filter_snr = SNR;
            
            [mx, idx2] = max(SNR); % the value with the most S/N
            
            obj.filter_loss = abs(mx-expected_snr)./mx;
%             obj.filter_loss = abs(mx-expected_snr);
            obj.filter_delta = sig(idx2)-sig(idx1);
            
            if idx2==1 || idx2==length(SNR)
                obj.filter_edge = 1;
            else
                obj.filter_edge = 0;
            end
            
        end
        
        function calcFit(obj, show_bit)
            
            import util.stat.sum2;
            
            obj.image_filt = filter2(obj.streak.psf, obj.image);
            
            if obj.use_trunc_mask
                obj.trunc_mask = obj.image_filt./sqrt(obj.streak.noise_var.*sum2(obj.streak.psf))>obj.streak.threshold;
            else
                obj.trunc_mask = false(size(obj.image_filt));
            end
            
            if obj.use_star_mask
                obj.star_mask = imdilate(obj.image_filt./sqrt(obj.streak.noise_var.*sum2(obj.streak.psf))>obj.streak.threshold*2, ones(9));
            else
                obj.trunc_mask = false(size(obj.image_filt));
            end
            
            % remove the mean (without stars)
            I = obj.image;
            I(obj.star_mask) = NaN;
            I = I - util.stat.median2(I);
            I(isnan(I)) = 0;
            obj.image = I;
            
            I_sigma_clipped = obj.image;
%             mask = obj.image>sqrt(obj.streak.noise_var);
            
%             neighborhood = true(3);
%             mask = imerode(mask, neighborhood);
%             mask = imdilate(mask, neighborhood);
%             I_sigma_clipped(~mask) = 0;
            
            obj.fitter = radon.Fitter;
            obj.fitter.show_bit = show_bit;
            obj.fitter.initial_x1 = obj.streak.x1_frt - obj.im_offset(2);
            obj.fitter.initial_x2 = obj.streak.x2_frt - obj.im_offset(2);
            obj.fitter.initial_y1 = obj.streak.y1_frt - obj.im_offset(1);
            obj.fitter.initial_y2 = obj.streak.y2_frt - obj.im_offset(1);
            obj.fitter.amplitude_normalization = obj.streak.I_frt;
            obj.fitter.psf_sigma = obj.streak.psf_sigma;
            obj.fitter.variance = obj.streak.noise_var;
            
            obj.fitter.run(I_sigma_clipped);
            
%             fprintf('FITTER: x1= %g | x2= %g | y1= %g | y2= %g\n', obj.fitter.x1, obj.fitter.x2, obj.fitter.y1, obj.fitter.y2);
%             fprintf('STREAK: x1= %g | x2= %g | y1= %g | y2= %g\n', obj.streak.x1, obj.streak.x2, obj.streak.y1, obj.streak.y2);
            
            obj.angle_correction = obj.streak.th - obj.fitter.angle;
            
            obj.streak.x1_fit = obj.fitter.x1 + obj.im_offset(2);
            obj.streak.x2_fit = obj.fitter.x2 + obj.im_offset(2);
            obj.streak.y1_fit = obj.fitter.y1 + obj.im_offset(1);
            obj.streak.y2_fit = obj.fitter.y2 + obj.im_offset(1);
            
            obj.streak.I_fit = obj.fitter.amplitude;
            obj.streak.snr_fit = obj.streak.I_fit./sqrt(2.*sqrt(pi).*obj.streak.psf_sigma.*obj.streak.noise_var./(obj.streak.L));
            
            obj.good_fitness = obj.fitter.chi2./obj.fitter.dof;
            
            if obj.checkFitGood % check the fit is even remotely successfull... 
                
                obj.streak.x1 = obj.streak.x1_fit;
                obj.streak.x2 = obj.streak.x2_fit;
                obj.streak.y1 = obj.streak.y1_fit;
                obj.streak.y2 = obj.streak.y2_fit;

                obj.streak.th = obj.fitter.angle;            
                obj.streak.L = sqrt((obj.streak.x2-obj.streak.x1).^2 + (obj.streak.y2-obj.streak.y1).^2);

                obj.streak.I = obj.streak.I_fit;
                obj.streak.snr = obj.streak.snr_fit;
            
            end
            
            
            
        end
        
        function rotateCrop(obj) % make a rotated image and masks
            
            import util.stat.sum2;
            
            obj.image_rot = imrotate(obj.image, obj.streak.th, 'bilinear','crop');
            norm = sqrt(obj.streak.noise_var.*sum2(obj.streak.psf));
            obj.trunc_mask_rot = filter2(obj.streak.psf, obj.image_rot)./norm>obj.streak.threshold;
            obj.star_mask_rot = imdilate(filter2(obj.streak.psf, obj.image_rot)./norm>obj.streak.threshold*2, ones(9));
            
            mid = floor(util.vec.imsize(obj.image)/2)+1;
            
%             obj.line_start = floor((obj.fitter.x1-mid(2)).*cosd(obj.fitter.angle) + (obj.fitter.y1-mid(1)).*sind(obj.fitter.angle) + mid(2));
            obj.line_start = floor((obj.streak.x1-obj.im_offset(2)-mid(2)).*cosd(obj.streak.th) + (obj.streak.y1-obj.im_offset(1)-mid(1)).*sind(obj.streak.th) + mid(2));
            
            if obj.line_start<1
                obj.line_start = 1;
            end
            
%             obj.line_end = ceil((obj.fitter.x2-mid(2)).*cosd(obj.fitter.angle) + (obj.fitter.y2-mid(1)).*sind(obj.fitter.angle) + mid(2));
            obj.line_end = ceil((obj.streak.x2-obj.im_offset(2)-mid(2)).*cosd(obj.streak.th) + (obj.streak.y2-obj.im_offset(1)-mid(1)).*sind(obj.streak.th) + mid(2));
            
            if obj.line_end<obj.line_start
                temp = obj.line_end;
                obj.line_end = obj.line_start;
                obj.line_start = temp;
            end
            
            if obj.line_end>size(obj.image,2)
                obj.line_end = size(obj.image,2);
            end
            
            % make a cropped image and masks
%             [~, pos] = max(sum(obj.image_rot,2));
            pos = round([obj.streak.midpoint_y, obj.streak.midpoint_x] - obj.im_offset);
            rel_pos = pos - mid;
            rel_pos = rel_pos*[cosd(obj.streak.th) sind(obj.streak.th); -sind(obj.streak.th) cosd(obj.streak.th)];
            pos = round(rel_pos + mid);
            
            top = pos(1) - floor(obj.width/2-1);
            bottom = pos(1) + floor(obj.width/2+1);
            
            if top<1 
                top = 1; 
                bottom = obj.width;
            elseif bottom>size(obj.image_rot,1) 
                top = size(obj.image_rot,1)-obj.width+1;
                bottom = size(obj.image_rot,1); 
            end
            
            obj.image_crop = obj.image_rot(top:bottom,:);
            obj.trunc_mask_crop = obj.trunc_mask_rot(top:bottom,:);
            obj.star_mask_crop = obj.star_mask_rot(top:bottom,:);
            
        end
        
        function val = checkFitGood(obj)
            
            val = true;
            
            val = val && (obj.fitter.chi2/obj.fitter.dof < 3);
            
            L = sqrt((obj.streak.x2_fit-obj.streak.x1_fit).^2+(obj.streak.y2_fit-obj.streak.y1_fit).^2);
            
            val = val && L<2*max(util.vec.imsize(obj.image)) && L>5;
            
            val = val && obj.fitter.chi2>0;
            
            val = val && obj.streak.snr_fit>0;
            
        end
        
        function calcBestAngle(obj) % not used
            
            import util.stat.sum2;
            
            angles = -obj.angle_range:obj.angle_step:obj.angle_range;
            
            score = zeros(length(angles),1);
            pos = zeros(length(angles),1);
            
            
            obj.image_filt = filter2(obj.streak.psf, obj.image);
            obj.trunc_mask = obj.image_filt./sqrt(obj.streak.noise_var.*sum2(obj.streak.psf))>obj.streak.threshold;
            obj.star_mask = imdilate(obj.image_filt./sqrt(obj.streak.noise_var.*sum2(obj.streak.psf))>obj.streak.threshold*2, ones(9));
            
            % remove the mean (without stars)
            I = obj.image;
            I(obj.star_mask) = NaN;
            I = I - util.stat.median2(I);
            I(isnan(I)) = 0;
            obj.image = I;
            
            I = obj.image_filt;
            I(obj.trunc_mask) = obj.streak.threshold.*sqrt(obj.streak.noise_var.*sum2(obj.streak.psf));
            I(obj.star_mask) = 0;
            
            for ii = 1:length(angles)
                
                Ir = imrotate(I, angles(ii)+obj.streak.th, 'bilinear','crop');
                [score(ii), pos(ii)] = max(sum(Ir(:,obj.l1:obj.l2),2));
                
                if obj.show_bit
                    plot(sum(Ir,2));
                    hold on;
                    xlabel(['angle= ' num2str(angles(ii)) ' | score= ' num2str(score(ii))]);
                    ax=gca;
                    ax.YLim = [0 1.5.*max(score(1))];
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
        
        function calcLength(obj) % not used
            
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
%             new_value = 2*pixel_width./obj.streak.psf_sigma;
            
            if pixel_width>obj.streak.subtract_psf_widths
                obj.streak.subtract_psf_widths = pixel_width;
            end
            
        end
        
        function calcEdgyness(obj)
            
            I = obj.image_crop;
            
            k = [1, 2, 1; 0, 0, 0; -1, -2, -1]; % sobel operator in Y
            
            If = filter2(k, I);
%             If2 = filter2(k',I);
%             If = sqrt(If.^2+If2.^2);
            
%             Sf = max(abs(sum(If(:,obj.line_start:obj.line_end),2)));
%             S = max(sum(I(:,obj.line_start:obj.line_end),2));

            obj.edge_norm = median(max(I(:,obj.line_start:obj.line_end), [], 1)); % average peak brightness
            
%             obj.edgyness = Sf./S;
%             obj.edgyness = obj.edgyness.*sqrt(obj.streak.psf_sigma);

            obj.edge_map = abs(If./obj.edge_norm);
            obj.edgyness = median(max(obj.edge_map(:, obj.line_start:obj.line_end),[],1));
            
%             obj.edgyness = obj.edgyness./20; % empirical normalization
            obj.edgyness = obj.edgyness.*sqrt(obj.streak.psf_sigma);
            
        end
        
        function calcVariability(obj)
            
            if obj.line_start+5>obj.line_end
                obj.variability = 0;
                obj.glint_power = 0;
                disp('streak found is too short for variability / glint calculation!');
                return;
            end
            
            a = obj.amplitudes(obj.line_start+2:obj.line_end-2);
            
            v = mean(a) + obj.streak.noise_var; % background and source noise together (assuming GAIN=1)
            
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
%             m1 = 0;
            m2 = sum(slice.*(y-m1).^2)./sum(slice);
            m3 = sum(slice.*(y-m1).^3)./sum(slice);
            
            obj.lateral_spread = sqrt(m2)./obj.streak.psf_sigma;
%             obj.lateral_skewness = m3./m2.^(3/2)./obj.streak.psf_sigma;
            obj.lateral_skewness = m3./obj.streak.psf_sigma;

            obj.lateral_profile = slice;
            
        end
        
        function calcBrightest(obj)
            
            obj.brightest_value = util.stat.max2(obj.image_crop);
            obj.brightest_number = nnz(obj.image_crop==obj.brightest_value);
            obj.brightest_ratio = obj.brightest_number./obj.streak.L;
            
        end
        
        function calcPeriodicity(obj) % not used
            
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
        
        function fitLine(obj) % not used
            
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
        
        function calcFiltering2(obj) % to be depricated
            
            import util.stat.sum2;
            import util.stat.max2;
            
            filter_range = -3.5:0.1:3.5;
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
%                 If(If./sqrt(obj.streak.noise_var)>obj.streak.threshold) = obj.streak.threshold.*sqrt(obj.streak.noise_var);
                If = util.img.truncate(If, obj.streak.threshold, obj.streak.noise_var);
                
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
        
        function calcFiltering3(obj) % to be depricated! 
            
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
        
        function correctEdges(obj) % to be depricated! 
            
            if abs(obj.angle_correction)<obj.angle_range
                obj.streak.th = obj.streak.th + obj.angle_correction; 
            end
            
            L = obj.streak.L;
            th = obj.streak.th;
            
            S = length(obj.amplitudes);
            mid_amp = floor(S/2)+1;
            
            distance_start = mid_amp - obj.line_start;
            distance_end = obj.line_end - mid_amp;
            
            delta_start = distance_start - L/2;
            delta_end = distance_end - L/2;
            
%             if obj.streak.th<0
%                 delta_start = -(L/2 - distance_end);
%                 delta_end = -(L/2 - distance_start);
%             end
            
            obj.streak.x1 = (obj.streak.x1 - cosd(th).*delta_start);
            obj.streak.y1 = (obj.streak.y1 - sind(th).*delta_start);
            
            if delta_start>0
                obj.streak.x1_subtraction = obj.streak.x1;
                obj.streak.y1_subtraction = obj.streak.y1;
            end
            
            obj.streak.x2 = (obj.streak.x2 + cosd(th).*delta_end);
            obj.streak.y2 = (obj.streak.y2 + sind(th).*delta_end);
            
            if delta_end>0
                obj.streak.x2_subtraction = obj.streak.x2;
                obj.streak.y2_subtraction = obj.streak.y2;
            end
            
            obj.streak.L = sqrt((obj.streak.x2-obj.streak.x1).^2 + (obj.streak.y2-obj.streak.y1).^2);
            
            
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
            line(obj.streak.psf_sigma_effective*[1 1], input.axes.YLim, 'LineStyle', '--', 'Color', 'green', 'Parent', input.axes);
            text(obj.streak.psf_sigma_effective*1.05, mean(input.axes.YLim), sprintf('\\sigma= %4.2f', obj.streak.psf_sigma_effective),...
                'Rotation', 90, 'VerticalAlignment', 'bottom', 'Color', 'green', 'Parent', input.axes);
            
            util.plot.inner_title(sprintf('loss= %5.3f', obj.filter_loss), 'Position', 'Bottom', 'axes', input.axes);
            
            if obj.filter_edge
                edge = '(out of bounds)';
            else
                edge = '';
            end
                
            util.plot.inner_title(sprintf('delta= %5.3f %s', obj.filter_delta, edge), 'Position', 'Top', 'axes', input.axes);
            
        end
        
        function plotAmplitudes(obj, varargin)
            
            if isempty(obj.amplitudes)
                return;
            end
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            prev_hold_state = input.axes.NextPlot;
            
            plot(input.axes, 1:length(obj.amplitudes), obj.amplitudes);
            
            if obj.line_start+5<=obj.line_end
            
                hold(input.axes, 'on');
                plot(input.axes, obj.line_start+2:obj.line_end-2, obj.amplitudes(obj.line_start+2:obj.line_end-2));
    %             plot(input.axes, obj.line_start:obj.line_end, lowpass(obj.amplitudes(obj.line_start:obj.line_end), obj.lowpass_freq));
                plot(input.axes, obj.line_start+2:obj.line_end-2, obj.amp_lowpass);
                input.axes.NextPlot = prev_hold_state;
                
            end
            
            xlabel(input.axes, 'position along streak');
            ylabel(input.axes, 'amplitude');
            
            util.plot.inner_title(['var= ' num2str(obj.variability)], 'Position', 'Bottom', 'axes', input.axes);
            util.plot.inner_title(['glint= ' num2str(obj.glint_power)], 'Position', 'Top', 'axes', input.axes);
            
        end
        
        function plotLateralProfile(obj, varargin)
            
            if isempty(obj.lateral_profile)
                return;
            end
            
            input = util.text.InputVars;            
            input.input_var('axes', [], 'axis');
            input.scan_vars(varargin);
            
            if isempty(input.axes)
                input.axes = gca;
            end
            
            x = (1:length(obj.lateral_profile))' - floor(length(obj.lateral_profile)/2) - 1;
            
            plot(input.axes, x, obj.lateral_profile);
            
            input.axes.XLim = [min(x), max(x)];
            
            xlabel(input.axes, 'lateral position');
            ylabel(input.axes, 'summed intensity');
            
            util.plot.inner_title(sprintf('spread= %4.2f | skew= %4.2f', obj.lateral_spread, abs(obj.lateral_skewness)), 'Position', 'Bottom', 'axes', input.axes);
            util.plot.inner_title(['edgyness= ' num2str(obj.edgyness)], 'Position', 'Top', 'axes', input.axes);
            
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

