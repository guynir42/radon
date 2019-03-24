classdef Classifier < handle
% This class takes streaks that had their Profiler activated, 
% and calculates scores on different metrics for classification as:
%   1) asteroids and comets (NEOs)
%   2) satellites in LEO or in high orbit
%   3) image artefacts and cosmic rays
%
% The metrics used for discrimination are:
% -Length of the streak.
% -filter size that gives maximum SNR. 
% -correction angle?
% -changing offsets (curved line).
% -large average spread.
% -periodic or changing amplitude. 
% -best fit low order polynomial?
% -

    properties(Transient=true)
        
    end
    
    properties % objects
        
    end
    
    properties % inputs/outputs
        
    end
    
    properties % switches/controls
        
        length_max = 200;
        filter_loss_max = 0.1;
        filter_delta_max = 0.5;
        edgyness_max = 5.0;
%         variability_max = 3; % not used
%         loss_var_product_max = 0.1; % not used 
        brightest_ratio_max = 0.5; 
        saturation_value = 6e4; % lower than 2^16 since we also allow for bias subtraction... 
        lateral_spread_max = 3.0;
        glint_power_max = 1.5;
        
        debug_bit = 1;
        
    end
    
    properties(Dependent=true)
        
        
        
    end
    
    properties(Hidden=true)
       
        version = 1.00;
        
    end
    
    methods % constructor
        
        function obj = Classifier(varargin)
            
            if ~isempty(varargin) && isa(varargin{1}, 'radon.Classifier')
                if obj.debug_bit, fprintf('Classifier copy-constructor v%4.2f\n', obj.version); end
                obj = util.oop.full_copy(varargin{1});
            else
                if obj.debug_bit, fprintf('Classifier constructor v%4.2f\n', obj.version); end
            
            end
            
        end
        
    end
    
    methods % reset/clear
        
    end
    
    methods % getters
        
    end
    
    methods % setters
        
    end
    
    methods % calculations
        
        function input(obj, streak, varargin)
            
            if streak.L>obj.length_max
                streak.type = 'too long';
                streak.addNote('Too long!');
                return;
            end
                
            if abs(mod(streak.th,90))<1
                streak.type = 'artefact';
                streak.addNote('Angle close to horizontal/vertical.');
                return;
            end
                
            if streak.prof.filter_edge
                streak.type = 'artefact';
                streak.addNote('Filter size is at edge of range.');
                return;
            end
            
            if streak.prof.filter_loss>obj.filter_loss_max 
                streak.type = 'artefact';
                streak.addNote(['Filter loss= ' num2str(streak.prof.filter_loss) '.']);
                return;
            end
            
            if streak.prof.filter_delta>obj.filter_delta_max
                streak.type = 'artefact';
                streak.addNote(['Filter delta= ' num2str(streak.prof.filter_delta) '.']);
                return;
            end
            
            if streak.prof.fitter.chi2==0
                streak.type = 'artefact';
                streak.addNote('Bad fit!');
                return;
            end
            
            if streak.prof.edgyness>obj.edgyness_max
                streak.type = 'artefact';
                streak.addNote(['Edgyness= ' num2str(streak.prof.edgyness) '.']);
                return;
            end
            
%             if streak.prof.variability>obj.variability_max
%                 streak.type = 'artefact';
%                 streak.addNote(['variablity= ' num2str(streak.prof.variability)]);
%                 return;
%             end
            
%             if streak.prof.variability.*streak.prof.filter_loss>obj.loss_var_product_max
%                 streak.type = 'artefact';
%                 streak.addNote(['variablity*loss= ' num2str(streak.prof.variability.*streak.prof.filter_loss)]);
%                 return;
%             end
            
            if streak.prof.brightest_ratio>obj.brightest_ratio_max
                if streak.prof.brightest_value>obj.saturation_value % saturated streak
                    streak.type = 'artefact';
                    streak.addNote(['More than ' num2str(100*streak.prof.brightest_ratio) '% of streak pixels are saturated!']);
                    return;
                else % this was probably close to the truncation threshold
                    streak.addNote(['Many points on this streak have the same maximal value: ' num2str(streak.prof.brightest_value) '.']);
                end
            end


            if streak.prof.lateral_spread>obj.lateral_spread_max % do brighter streaks also have higher spreads??
                streak.type = 'artefact';
                streak.addNote(['Lateral spread= ' num2str(streak.prof.lateral_spread) '.']);
                return;
            end
            
            if streak.prof.glint_power>obj.glint_power_max
                streak.type = 'satellite';
                streak.addNote(['Glint= ' num2str(streak.prof.glint_power) '.']);
                return;
            end
            
            % if all other tests do not fail, we must assume it is an asteroid
            streak.type = 'asteroid';
            
        end
        
    end
    
    methods % plotting tools / GUI
        
        function printStats(obj, streaks)
            
            for ii = 1:length(streaks)
                
                f = streaks(ii).frame_num;
                type = streaks(ii).type;
                SNR = streaks(ii).snr;
                L = round(streaks(ii).L);
                loss = streaks(ii).prof.filter_loss;
                edgy = streaks(ii).prof.edgyness;
                star = ' '; 
                if streaks(ii).prof.filter_edge, star = '*'; end
                vari = streaks(ii).prof.variability;
                glint = streaks(ii).prof.glint_power;
                spre = streaks(ii).prof.lateral_spread;
                
%                 fprintf('ii= % 3d f= % 3d (%-9s) | S/N= %6.2f | L= % 4d | loss= %6.4f%s | edgy= %6.4f | var= %6.4f | loss*var*10= %4.2f | glint= %4.2f | spread= %4.2f \n',...
%                     ii, f, type, SNR, L, loss, star, edgy, vari, loss*vari*10, glint, spre);
                fprintf('ii= % 3d f= % 3d (%-9s) | S/N= %6.2f | L= % 4d | loss= %6.4f%s | edgy= %6.4f | glint= %4.2f | spread= %4.2f \n',...
                    ii, f, type, SNR, L, loss, star, edgy, glint, spre);
                
            end
            
        end
        
    end    
    
end

