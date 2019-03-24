classdef Fitter < handle
    
    properties % fit results
        
        bias;
        amplitude;        
        x1;
        x2;
        y1;
        y2;
        sigma;
        
        chi2;
        dof;
        p_value;
        
        model; % last calculated model image
        
    end
    
    properties % starting parameters / inputs
    
        image;
        variance = 1; % map or scalar
        gain = 1; % this helps figure out if the fit is good
        psf_sigma = 1;
        
        initial_x1;
        initial_x2;
        initial_y1;
        initial_y2;
        
    end
    
    properties % fit options
        
        use_par_mapping = 0;
        use_photometry = 1;
        
        mode = 'minimize'; % can choose "minimize" or "scan"
        use_sigma_fitting = 0;
        use_priors = 1;
        
        range_amplitude = 1;
        range_angle = 20;
        range_midpoint = 15;
        range_sigma = 1;
        
%         options = {'TolFun', 1, 'TolX', 0.1};
        options = {};
        
        show_bit = 0;

    end
        
    properties (Hidden=true)
        
        amplitude_normalization;
        
        multiply_range = 3;
        penalty_power = 30;
        
    end
    
    methods % calculations
        
        function obj = Fitter(obj, varargin) % constructor
            
            % allow constructor to parse inputs like initial parameters and
            % fitting options. 
            
            % ...
            
        end
        
        function clear(obj)

            obj.image = [];
            
            obj.amplitude = [];
            obj.sigma = [];
            obj.x1 = [];
            obj.x2 = [];
            obj.y1 = [];
            obj.y2 = [];

            obj.chi2 = [];
            
        end
        
        function run(obj, I)
            
            if nargin>1 && ~isempty(I)
                obj.image = I;
            end
            
            if util.text.cs(obj.mode, 'minimizer')
                obj.runMinimizer;
            elseif util.text.cs(obj.mode, 'scan')
                obj.runScan;
            else
                error('Unknown fitting option "%s". Use "minimizer" or "scan"...', obj.mode);
            end
            
            if obj.use_photometry
                obj.photometry;
            end
            
            obj.getChi2;
            
        end
       
        function runMinimizer(obj)
            
            import util.stat.sum2;
                
            im_size = util.vec.imsize(obj.image);
            
            % normalize initial values to something close to unity
            bias_norm = util.stat.corner_median(obj.image);
%             L = sqrt((obj.initial_x2-obj.initial_x1).^2 + (obj.initial_y2-obj.initial_y1).^2);
%             obj.amplitude_normalization = sum2(obj.image)/L;
            amp_norm = 1; % normalize the initial amplitude to 1
            
            if obj.use_par_mapping % use more parameters that are less correlated
                
                mid_x = mean([obj.initial_x1,obj.initial_x2]);
                mid_y = mean([obj.initial_y1,obj.initial_y2]);
                r1 = sqrt((obj.initial_x1-mid_x).^2+(obj.initial_y1-mid_y).^2);
                r2 = sqrt((obj.initial_x2-mid_x).^2+(obj.initial_y2-mid_y).^2);
                
                mid_x_norm = mid_x./im_size(2);
                mid_y_norm = mid_y./im_size(1);
                r1_norm = r1./sqrt(prod(im_size));
                r2_norm = r2./sqrt(prod(im_size));
                
                angle = atan2d(obj.initial_y2-obj.initial_y1,obj.initial_x2-obj.initial_x1);
                angle_norm = angle./90;
                
                b_initial = [mid_x_norm, mid_y_norm, r1_norm, r2_norm, angle_norm, amp_norm, bias_norm];
                
            else
                
                x1_norm = obj.initial_x1./im_size(2); % normalize to relative units
                x2_norm = obj.initial_x2./im_size(2); % normalize to relative units
                y1_norm = obj.initial_y1./im_size(1); % normalize to relative units
                y2_norm = obj.initial_y2./im_size(1); % normalize to relative units

                b_initial = [x1_norm, x2_norm, y1_norm, y2_norm, amp_norm, bias_norm];
            
            end
            
            if iscell(obj.options)
                opt = optimset(obj.options{:});
            else
                opt = obj.options;
            end

            func = @(b) obj.calcMatch(b);

            if obj.use_sigma_fitting
                b_initial = [b_initial, obj.psf_sigma];
            end

            b_found = fminsearch(func, b_initial, opt);

            % translate back to coordinates
            if obj.use_par_mapping
                mid_x = b_found(1).*im_size(2);
                mid_y = b_found(2).*im_size(1);
                r1 = b_found(3).*sqrt(prod(im_size));
                r2 = b_found(4).*sqrt(prod(im_size));
                angle = b_found(5).*90;
                
                obj.x1 = mid_x-r1.*cosd(angle);
                obj.y1 = mid_y-r1.*sind(angle);
                obj.x2 = mid_x+r2.*cosd(angle);
                obj.y2 = mid_y+r2.*sind(angle);
                
                obj.amplitude = b_found(6).*obj.amplitude_normalization; 
                obj.bias = b_found(7); % no normalization
                
            else
                
                obj.x1 = b_found(1).*im_size(2);
                obj.x2 = b_found(2).*im_size(2);
                obj.y1 = b_found(3).*im_size(1);
                obj.y2 = b_found(4).*im_size(1);
                
                obj.amplitude = b_found(5).*obj.amplitude_normalization; 
                obj.bias = b_found(6); % no normalization
                
            end
            
            if obj.use_sigma_fitting
                obj.sigma = b_found(end);
            else
                obj.sigma = obj.psf_sigma;
            end
            
        end
        
        function [chi2, dof] = calcMatch(obj, b_vec_norm)
            
            im_size = util.vec.imsize(obj.image);
            
            if obj.use_par_mapping
                
                mid_x = b_vec_norm(1).*im_size(2);
                mid_y = b_vec_norm(2).*im_size(1);
                r1 = b_vec_norm(3).*sqrt(prod(im_size));
                r2 = b_vec_norm(4).*sqrt(prod(im_size));
                angle = b_vec_norm(5).*90;
                amp = b_vec_norm(6).*obj.amplitude_normalization; 
                bias = b_vec_norm(7);
                
                x1 = mid_x-r1.*cosd(angle);
                y1 = mid_y-r1.*sind(angle);
                x2 = mid_x+r2.*cosd(angle);
                y2 = mid_y+r2.*sind(angle);
                
            else
                
                x1 = b_vec_norm(1).*im_size(2);
                x2 = b_vec_norm(2).*im_size(2);
                y1 = b_vec_norm(3).*im_size(1);
                y2 = b_vec_norm(4).*im_size(1); 
                amp = b_vec_norm(5).*obj.amplitude_normalization; 
                bias = b_vec_norm(6); 
                angle = atan2d(y2-y1, x2-x1);
                mid_x = mean([x1,x2]);
                mid_y = mean([y1,y2]);
                
            end
            
            if obj.use_sigma_fitting
                sigma = b_vec_norm(end); % sigma
            else
                sigma = obj.psf_sigma;
            end
            
            % calculate prior penalty for going too far off initial values
            prior = 1;
            
            % amplitude penalty...
%             prior = prior.*cosh((b_vec_norm(5)-1)./3./obj.multiply_range).^obj.penalty_power;
%             prior = prior.*(1+b_vec_norm(5).^2)./b_vec_norm(5).^2; % make sure amplitude doesn't go to zero...
%             disp(['b_vec_norm(5)= ' num2str(b_vec_norm(5)) ' | prior= ' num2str(prior)]);
            
            % angle penalties...
%             initial_angle = atan2d(obj.initial_y2-obj.initial_y1, obj.initial_x2-obj.initial_x1); 
%             prior = prior.*cosh((angle-initial_angle)./obj.range_angle./obj.multiply_range).^obj.penalty_power;
            
            % distance from line penalty...
            if obj.use_par_mapping
            
                initial_mid_x = mean([obj.initial_x1,obj.initial_x2]);
                initial_mid_y = mean([obj.initial_y1,obj.initial_y2]);
                dist = sqrt((mid_x-initial_mid_x).^2 + (mid_y-initial_mid_y).^2);
                prior = prior.*cosh(dist./obj.range_midpoint./obj.multiply_range).^obj.penalty_power;
                
            end

%             relative_angle = atan2d(mid_y-initial_mid_y,mid_x-initial_mid_x)-initial_angle; % between midpoints, relative to initial position
%             midpoint_distance = sqrt((initial_mid_y-mid_y).^2 + (initial_mid_x-mid_x).^2);
%             relative_distance = midpoint_distance.*sind(relative_angle);
%             prior = prior.*cosh(relative_distance./obj.range_midpoint./obj.multiply_range).^obj.penalty_power;

            model = bias + amp.*radon.model(im_size, x1, x2, y1, y2, sigma, [], [], 0); % last zero is for no downsampling! 
            obj.model = model;
            diff_image = obj.image - model;
            
            chi2 = util.stat.sum2( (diff_image).^2./obj.variance );
            
            if obj.use_priors
                chi2 = chi2.*prior;
            end
            
            dof = nnz(~isnan(diff_image)) - length(b_vec_norm);
            
            if obj.show_bit
                util.plot.show(diff_image, 'autodyn', 1);
                xlabel(['chi2= ' num2str(chi2) ' | dof= ' num2str(dof)]);
                title(num2str(b_vec_norm));
                
%                 fprintf('x1= %f x2= %f y1= %f y2= %f\n', x1,x2,y1,y2);
                drawnow;
            end
            
        end
        
        function [chi2, dof] = calcMatchOld(obj, b_vec_norm) % to be depricated!
            
            im_size = util.vec.imsize(obj.image);
            
            b_vec(1) = b_vec_norm(1); 
            b_vec(2) = b_vec_norm(2).*obj.amplitude_normalization;
            b_vec(3) = b_vec_norm(3).*im_size(2);
            b_vec(4) = b_vec_norm(4).*im_size(1);
            b_vec(5) = b_vec_norm(5).*90;
            b_vec(6) = b_vec_norm(6).*sqrt(prod(im_size));
            b_vec(7) = b_vec_norm(7).*sqrt(prod(im_size));
            if length(b_vec_norm)>=8
                b_vec(8) = b_vec_norm(8);
            end
            
            % calculate prior penalty for going to far
            prior = 1;
            
            % amplitude penalty...
%             prior = prior.*cosh((b_vec_norm(2)-1)./3./obj.multiply_range).^obj.penalty_power;
%             prior = prior.*(1+b_vec_norm(2).^2)./b_vec_norm(2).^2; % make sure amplitude doesn't go to zero...
%             disp(['b_vec_norm(2)= ' num2str(b_vec_norm(2)) ' | prior= ' num2str(prior)]);
            
            % midpoint penalties...
            initial_mid_x = mean([obj.initial_x1,obj.initial_x2]);
            initial_mid_y = mean([obj.initial_y1,obj.initial_y2]);
            prior = prior.*cosh((b_vec(3)-initial_mid_x)./obj.range_midpoint./obj.multiply_range).^obj.penalty_power;
            prior = prior.*cosh((b_vec(4)-initial_mid_y)./obj.range_midpoint./obj.multiply_range).^obj.penalty_power;
            
            % angle penalties...
            initial_angle = atan2d(obj.initial_y2-obj.initial_y1, obj.initial_x2-obj.initial_x1);
            prior = prior.*cosh((b_vec(5)-initial_angle)./obj.range_angle./obj.multiply_range).^obj.penalty_power;
            
            if obj.use_sigma_fitting
                sigma = b_vec(8);
            else
                sigma = obj.psf_sigma;
            end
            
            model = radon.Fitter.getModel(im_size, b_vec(3), b_vec(4), b_vec(5), b_vec(6), b_vec(7), sigma);
            model = b_vec(1) + b_vec(2).*model;
            diff_image = obj.image - model;
            
            chi2 = util.stat.sum2( (diff_image).^2./obj.variance );
            
            if obj.use_priors
                chi2 = chi2.*prior;
            end
            
            dof = nnz(~isnan(diff_image)) - length(b_vec);
            
            if obj.show_bit
                util.plot.show(diff_image);
                xlabel(['chi2= ' num2str(chi2) ' | dof= ' num2str(dof)]);
                title(num2str(b_vec));
                drawnow;
            end
            
        end
        
        function runScan(obj)
            
        end
        
        function val = angle(obj)
            
            val = atan2d(obj.y2-obj.y1, obj.x2-obj.x1);
            
        end
        
        function photometry(obj) % fit only amplitude and angle, with NaNs far away from the streak. Also goodness of fit. 
            
            im_size = util.vec.imsize(obj.image);
            
            model = @(b) b(1).*obj.amplitude.*radon.model(im_size, obj.x1, obj.x2, obj.y1, obj.y2, obj.sigma, NaN, [], 2) + b(2);
            
            diff = @(b) obj.image - model(b);
            
            func = @(b) util.stat.sum2(diff(b).^2); % use the default threshold and oversample
            
            if iscell(obj.options)
                opt = optimset(obj.options{:});
            else
                opt = obj.options;
            end
            
            b_found = fminsearch(func, [1, obj.bias], opt); 
            
%             M = model(b_found);
            
%             bg_var = util.stat.sum2(~isnan(M).*obj.variance); % number of pixels used in the fit times their variance
%             L = sqrt((obj.x2-obj.x1).^2+(obj.y2-obj.y1).^2);
%             source_var = L*b_found(1); % source noise var is proportional to the amplitude 
            
%             total_std = sqrt(bg_var+source_var);
%             obj.norm_residuals = sqrt(util.stat.sum2((M-obj.image).^2))./total_std;
            
            % add test before updating amplitude? 
            obj.amplitude = b_found(1).*obj.amplitude; 
            obj.bias = b_found(2);
            
            
        end
        
        function [chi2, dof] = getChi2(obj)
            
            M = obj.getModel(NaN);
            
            obj.dof = nnz(~isnan(M)) - 6 - obj.use_sigma_fitting; 
%             bg_var = util.stat.sum2(~isnan(M).*obj.variance); % number of pixels used in the fit times their variance
%             L = sqrt((obj.x2-obj.x1).^2+(obj.y2-obj.y1).^2);
%             source_var = L*obj.amplitude; % source noise var is proportional to the amplitude 
            
            total_var = obj.variance+obj.gain.*M;
            obj.chi2 = util.stat.sum2((M-obj.image).^2./total_var);
            obj.p_value = chi2cdf(obj.chi2,obj.dof,'upper');
            
            if nargout>0
                chi2 = obj.chi2;
                dof = obj.dof;
                p_value = obj.p_value;
            end
            
        end
        
        function M = getModel(obj, replace_value, threshold, oversample)
            
            if nargin<2 || isempty(replace_value)
                replace_value = [];
            end
            
            if nargin<3 || isempty(threshold)
                threshold = [];
            end
            
            if nargin<4 || isempty(oversample)
                oversample = [];
            end
            
            M = obj.bias + obj.amplitude.*radon.model(util.vec.imsize(obj.image), obj.x1, obj.x2, obj.y1, obj.y2, obj.sigma, ...
                replace_value, threshold, oversample);
            
        end
        
    end

    methods % plotting tools
        
        
        
    end
    
    methods (Static=true)
        
        function M = getModelOld(im_size, mid_x, mid_y, angle, r1, r2, sigma)
            
            x1 = mid_x - r1.*cosd(angle);
            x2 = mid_x + r2.*cosd(angle);
            y1 = mid_y - r1.*sind(angle);
            y2 = mid_y + r2.*sind(angle);
            
            M = radon.model(im_size, x1, x2, y1, y2, sigma, 0, [], 0);
            
%             S = util.stat.sum2(M);
%             if S>0
%                 M = M./S;
%             end
            
        end
        
    end
    
end






