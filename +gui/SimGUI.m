classdef SimGUI < handle
    
    properties 
        
        owner@radon.Simulator; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
                
        panel_stats;
        input_im_size;
        input_intensity;
        input_bg_noise_var;
        button_use_source_noise;
        
        panel_psf;
        button_use_psf;
        input_psf_sigma;
        
        panel_random;
        button_use_random;
        input_num_streaks;
        input_intensity_range;
        input_length_range;
        input_midpoint_x_range;
        input_midpoint_y_range;
        input_theta_range;
        
        panel_control;
        choose_display_what;
        button_pick_points;
        button_run;
        button_finder;
        button_use_finder;
        
        panel_close;
        button_close;
        
        panel_image;
%         button_reset_axes;
        axes_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.01;
        
    end
            
    methods % constructor
       
        function obj = SimGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'radon.Simulator')
                
                if obj.debug_bit, fprintf('SimGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            elseif isa(owner, 'radon.gui.SimGUI')
                
                if obj.debug_bit, fprintf('SimGUI copy-constructor v%4.2f\n', obj.version); end
                
                obj = util.oop.full_copy(owner);
                
            else
                error('Input an radon.Simulator to constructor of SimGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function makeGUI(obj)
            
            import util.plot.GraphicButton;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('streak simulator');
            end
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            obj.fig.name = 'streak simulator';
            movegui(obj.fig.fig, 'center');
            total_num = 15;
                                   
            %%%%%%%%%%% panel stats %%%%%%%%%%%%%%%%%%
            
            num = 4; obj.panel_stats = uipanel('Title', 'statistics', 'Position', [0 (total_num-num)/total_num 0.2 num/total_num]);
            obj.input_intensity = GraphicButton(obj.panel_stats, [0 (num-1)/num 1 1/num], obj.owner, 'intensity', 'input', 'I= ');
            obj.input_im_size = GraphicButton(obj.panel_stats, [0 (num-2)/num 1 1/num], obj.owner, 'im_size', 'input', 'size= ');
            obj.input_bg_noise_var = GraphicButton(obj.panel_stats, [0 (num-3)/num 1 1/num], obj.owner, 'bg_noise_var', 'input', 'bg var= ');
            obj.button_use_source_noise = GraphicButton(obj.panel_stats, [0 (num-4)/num 1 1/num], obj.owner, 'use_source_noise', 'toggle', 'no source noise', 'source noise');
                        
            %%%%%%%%%%% panel PSF %%%%%%%%%%%%%%%%%%%%
            
            num = 2; obj.panel_psf = uipanel('Title', 'PSF', 'Position', [0 (total_num-num-4)/total_num 0.2 num/total_num]);
            obj.button_use_psf = GraphicButton(obj.panel_psf, [0 (num-1)/num 1 1/num], obj.owner, 'use_psf', 'toggle', 'no PSF', 'use PSF');
            obj.input_psf_sigma = GraphicButton(obj.panel_psf, [0 (num-2)/num 1 1/num], obj.owner, 'psf_sigma', 'input', 'width= ');
            
            %%%%%%%%%%% panel random %%%%%%%%%%%%%%%%%
            
            num = 7; obj.panel_random = uipanel('Title', 'randomize', 'Position', [0 (total_num-num-6)/total_num 0.2 num/total_num]);
            obj.button_use_random = GraphicButton(obj.panel_random, [0 (num-1)/num 1 1/num], obj.owner, 'use_random', 'toggle', 'no random', 'use random');
            obj.input_num_streaks = GraphicButton(obj.panel_random, [0 (num-2)/num 1 1/num], obj.owner, 'num_random_streaks', 'input', 'N= ');
            obj.input_intensity_range = GraphicButton(obj.panel_random, [0 (num-3)/num 1 1/num], obj.owner, 'intensity_range', 'input', 'I= ');
            obj.input_length_range = GraphicButton(obj.panel_random, [0 (num-4)/num 1 1/num], obj.owner, 'length_range', 'input', 'length= ');
            obj.input_midpoint_x_range = GraphicButton(obj.panel_random, [0 (num-5)/num 1 1/num], obj.owner, 'midpoint_x_range', 'input', 'x range= ');
            obj.input_midpoint_y_range = GraphicButton(obj.panel_random, [0 (num-6)/num 1 1/num], obj.owner, 'midpoint_y_range', 'input', 'y range= ');
            obj.input_theta_range = GraphicButton(obj.panel_random, [0 (num-7)/num 1 1/num], obj.owner, 'angle_range', 'input', 'angle= ');
            
            %%%%%%%%%%% panel control %%%%%%%%%%%%%%%%
            
            num = 6; obj.panel_control = uipanel('Title', 'control', 'Position', [0.2 0 0.8 2/total_num]);
            obj.choose_display_what = uicontrol(obj.panel_control, 'Style', 'popupmenu', 'Units', 'Normalized', 'Position', [0/num 0 1/num 1],...
                'Callback', @obj.callback_display_what, 'String', obj.owner.display_what_list);
            obj.button_pick_points = GraphicButton(obj.panel_control, [1/num 0 1/num 1], obj.owner, 'pickPoints', 'push', 'pick');
            obj.button_run = GraphicButton(obj.panel_control, [2/num 0 2.5/num 1], obj.owner, 'run', 'push', 'RUN');
            obj.button_finder = GraphicButton(obj.panel_control, [4.5/num 0 1/num 1], obj.owner, 'finder', 'push', 'Finder');
            obj.button_use_finder = GraphicButton(obj.panel_control, [5.5/num 0 0.5/num 1], obj.owner, 'use_finder', 'toggle', 'off', 'on');            
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Title', 'Close', 'Position', [0 0 0.2 2/total_num]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
                      
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
              
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 2/total_num 0.8 (total_num-2)/total_num]);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.update;
            
        end
        
        function makeAxes(obj, ~, ~)
            
            delete(obj.axes_image);
            
            obj.axes_image = axes('Parent', obj.panel_image);
            
            obj.panel_contrast.ax = obj.axes_image;
%             colorbar(obj.axes_image);
            axis(obj.axes_image, 'image');
            
            obj.panel_contrast.ax = obj.axes_image;
            
        end
                
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.choose_display_what.FontSize = obj.font_size;
            
        end
        
        function updateGUI(obj)
            obj.update;
        end
                
        function c = check(obj)
           
            c = ~isempty(obj.panel_image) && isvalid(obj.panel_image);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_display_what(obj, hndl, ~)
            
%             obj.owner.display_what_index = hndl.Value;
            obj.owner.display_what = obj.owner.display_what_list{hndl.Value};
            obj.owner.show;
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end