classdef FinderGUI < handle
    
    properties 
        
        owner@radon.Finder; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_prep;
        button_use_prep;
        button_use_sub;
        button_use_mean;
        button_use_median;
        button_use_conv;
        button_use_crop;
        input_crop_size;
        
        panel_search;
        input_thresh;
        button_use_short;
        input_min_length;
        button_use_exclude;
        input_exclude_dx;
        input_exclude_dy;
        
        button_use_save;
        
        panel_info;
        button_filename;
        button_batch;
        button_frame;
        
        panel_display;        
        button_close;
        input_line_offset;
        input_rect_size;
        button_use_monochrome;
        button_display_which;
        button_play;
        button_prev;
        input_streak_num;
        button_next;
        button_best;
        button_show;
        
        panel_image;
    
    end
    
    properties (Hidden=true)
              
        version = 1.02;
        
    end
            
    methods % constructor
       
        function obj = FinderGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'radon.Finder')
                
                if obj.debug_bit, fprintf('FinderGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            elseif isa(owner, 'radon.gui.FinderGUI')
                
                
                
            else
                error('Input an radon.Finder to constructor of FinderGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function makeGUI(obj)
            
            import util.plot.GraphicButton;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            if isempty(obj.fig)
                obj.fig = util.plot.FigHandler('FinderGUI');
            end
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            obj.fig.name = 'FinderGUI';
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel prep %%%%%%%%%%%%%%%%
            
            obj.panel_prep = uipanel('Title', 'preprocess', 'Position', [0.0 8/12 0.2 4/12]);
            
            num = 4;
                        
            obj.button_use_prep = GraphicButton(obj.panel_prep, [0 (num-1)/num 1 1/num], obj.owner, 'use_preprocess', 'toggle', 'no prep', 'use prep','','','','red');
            obj.button_use_mean = GraphicButton(obj.panel_prep, [0 (num-2)/num 0.5 1/num], obj.owner, 'use_subtract_mean', 'toggle', 'sub-mean', 'sub-mean','','','','red');            
            obj.button_use_median = GraphicButton(obj.panel_prep, [0.5 (num-2)/num 0.5 1/num], obj.owner, 'use_subtract_median', 'toggle', 'sub-median', 'sub-median','','','','red');
            obj.button_use_crop = GraphicButton(obj.panel_prep, [0 (num-3)/num 1 1/num], obj.owner, 'use_crop', 'toggle', 'no crop', 'use crop','','','','red');
            obj.input_crop_size = GraphicButton(obj.panel_prep, [0 (num-4)/num 1 1/num], obj.owner, 'crop_size', 'input', 'crop= ');
            
            %%%%%%%%%%% panel search %%%%%%%%%%%%%%%%
            
            obj.panel_search = uipanel('Title', 'search', 'Position', [0.0 2/12 0.2 6/12]);
                        
            num = 6;
            
            obj.input_thresh = GraphicButton(obj.panel_search, [0 (num-1)/num 1 1/num], obj.owner, 'threshold', 'input', 'threshold= ');
            obj.button_use_short = GraphicButton(obj.panel_search, [0 (num-2)/num 1 1/num], obj.owner, 'use_short', 'toggle', 'no short', 'use short','','','','red');
            obj.input_min_length = GraphicButton(obj.panel_search, [0 (num-3)/num 1 1/num], obj.owner, 'min_length', 'input', 'len= ');
            obj.button_use_exclude = GraphicButton(obj.panel_search, [0 (num-4)/num 1 1/num], obj.owner, 'use_exclude', 'toggle', 'no exclude', 'use exclude','','','','red');
            obj.input_exclude_dx = GraphicButton(obj.panel_search, [0 (num-5)/num 1 1/num], obj.owner, 'exclude_dx', 'input', 'exc. dx= ');
            obj.input_exclude_dy = GraphicButton(obj.panel_search, [0 (num-6)/num 1 1/num], obj.owner, 'exclude_dy', 'input', 'exc. dy= ');
            
%             obj.button_use_save = GraphicButton(obj.panel_prep, [0 (num-11)/num 1 1/num], obj.owner, 'use_save_images', 'toggle', 'no save', 'save images');
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = uipanel('Title', 'info', 'Position', [0.2 11/12 0.8 1/12]);
            
            obj.button_filename = GraphicButton(obj.panel_info, [0.0 0 0.8 1], obj.owner, 'filename', 'info');
%             obj.button_filename.control.HorizontalAlignment = 'right';
            obj.button_filename.font_size = 'small';
            obj.button_filename.Callback = @obj.callback_filename;
            obj.button_batch = GraphicButton(obj.panel_info, [0.8 0 0.1 1], obj.owner, 'batch_num', 'info', 'b= ');
            obj.button_frame = GraphicButton(obj.panel_info, [0.9 0 0.1 1], obj.owner, 'frame_num', 'info', 'f= ');
            
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            obj.panel_display = uipanel('Title', 'display', 'Position', [0.0 0.0 1 2/12]);
            
            obj.button_close = GraphicButton(obj.panel_display, [0.0 0.0 0.2 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.input_line_offset = GraphicButton(obj.panel_display, [0.2 0.5 0.2 0.5], obj.owner, 'line_offset', 'input', ' line offset= ');
            obj.input_rect_size = GraphicButton(obj.panel_display, [0.2 0.0 0.2 0.5], obj.owner, 'rect_size', 'input', ' rect size= ');
            
            obj.button_use_monochrome = GraphicButton(obj.panel_display, [0.4 0.5 0.2 0.5], obj.owner, 'display_monochrome', 'toggle', 'color', 'monochrome');
            
            obj.button_display_which = GraphicButton(obj.panel_display, [0.6 0.5 0.2 0.5], obj.owner, 'cycleDisplayWhich', 'push');
            obj.button_play = GraphicButton(obj.panel_display, [0.4 0 0.2 0.5], obj.owner, 'showAllStreaks', 'push', 'play');
            obj.button_prev = GraphicButton(obj.panel_display, [0.6 0 0.05 0.5], obj.owner, 'prevDisplayIndex', 'push', '-');
            obj.input_streak_num = GraphicButton(obj.panel_display, [0.65 0 0.05 0.5], obj.owner, 'display_index', 'input', ' ');
            obj.button_next = GraphicButton(obj.panel_display, [0.70 0 0.05 0.5], obj.owner, 'nextDisplayIndex', 'push', '+');
            obj.button_best = GraphicButton(obj.panel_display, [0.75 0 0.05 0.5], obj.owner, 'resetDisplayIndex', 'push', 'best');
            
            obj.button_show = GraphicButton(obj.panel_display, [0.8 0.0 0.2 1], obj.owner, 'show', 'push');
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 2/12 0.8 9/12]);
            
            obj.update;
            
        end
                            
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.button_display_which.String = obj.owner.display_which;
            
            obj.input_streak_num.String = obj.owner.getStreakIndex;
            
        end
            
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_image) && isvalid(obj.panel_image);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_filename(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: filename'); end
            
            disp(['filename: ' obj.owner.filename]);
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end