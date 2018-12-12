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
        filter_loss_max = 0.2;
        variability_max = 3;
        loss_var_product_max = 0.1;
        glint_power_max = 2.5;
        lateral_spread_max = 2;
        
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
                streak.type = 'satellite';
                streak.addNote('too long!');
                return;
            end
            
            if streak.prof.filter_edge
                streak.type = 'artefact';
                streak.addNote('filter size is at edge of range');
                return;
            end
            
            if streak.prof.filter_loss>obj.filter_loss_max 
                streak.type = 'artefact';
                streak.addNote(['filter loss= ' num2str(streak.prof.filter_loss)]);
                return;
            end
            
            if streak.prof.variability>obj.variability_max
                streak.type = 'artefact';
                streak.addNote(['variablity= ' num2str(streak.prof.variability)]);
                return;
            end
            
            if streak.prof.variability.*streak.prof.filter_loss>obj.loss_var_product_max
                streak.type = 'artefact';
                streak.addNote(['variablity*loss= ' num2str(streak.prof.variability.*streak.prof.filter_loss)]);
                return;
            end
                
            if streak.prof.glint_power>obj.glint_power_max
                streak.type = 'artefact';
                streak.addNote(['glint= ' num2str(streak.prof.glint_power)]);
                return;
            end
                        
            if streak.prof.lateral_spread>obj.lateral_spread_max % do brighter streaks also have higher spreads??
                streak.type = 'artefact';
                streak.addNote(['lateral spread= ' num2str(streak.prof.lateral_spread)]);
                return;
            end
            
            % add tests for amplitudes, offsets, spreads
            
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
                star = ' '; 
                if streaks(ii).prof.filter_edge, star = '*'; end
                vari = streaks(ii).prof.variability;
                glint = streaks(ii).prof.glint_power;
                spre = streaks(ii).prof.lateral_spread;
                
                fprintf('ii= % 3d f= % 3d (%-9s) | S/N= %6.2f | L= % 4d | loss= %6.4f%s | var= %6.4f | loss*var*10= %4.2f | glint= %4.2f | spread= %4.2f \n', ii, f, type, SNR, L, loss, star, vari, loss*vari*10, glint, spre);
                
            end
            
        end
        
    end    
    
end

