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
        
        panel_options;
        
        button_use_sub;
        
        button_use_autothresh;
        button_autothresh_val;
        input_thresh;
        
        button_use_short;
        input_min_length;
        input_pixel_res;
        
        button_use_recursive;
        input_recursion_depth;
        
        panel_exclude;
        button_use_exclude;
        input_exclude_dx;
        input_exclude_dy;
        
        button_use_crop;
        input_crop_size;
        
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
            
            %%%%%%%%%%% panel options %%%%%%%%%%%%%%%%
            
            obj.panel_options = uipanel('Title', 'options', 'Position', [0.0 0.2 0.2 0.8]);
            
            num = 14;
            
            obj.button_use_sub = GraphicButton(obj.panel_options, [0 (num-1)/num 1 1/num], obj.owner, 'use_subtract_mean', 'toggle', 'no sub. mean', 'use sub. mean');
            obj.button_use_autothresh = GraphicButton(obj.panel_options, [0 (num-2)/num 0.8 1/num], obj.owner, 'use_autothresh', 'toggle', 'no autothresh', 'use autothresh');
            obj.button_autothresh_val = GraphicButton(obj.panel_options, [0.8 (num-2)/num 0.2 1/num], obj.owner, 'latest_thresh', 'info', '');
            obj.button_autothresh_val.font_size = 'small';
            obj.input_thresh = GraphicButton(obj.panel_options, [0 (num-3)/num 1 1/num], obj.owner, 'threshold', 'input', 'threshold= ');
            
            obj.button_use_short = GraphicButton(obj.panel_options, [0 (num-4)/num 1 1/num], obj.owner, 'use_short', 'toggle', 'no short', 'use short');
            obj.input_min_length = GraphicButton(obj.panel_options, [0 (num-5)/num 1 1/num], obj.owner, 'min_length', 'input', 'len= ');
           
            obj.button_use_recursive = GraphicButton(obj.panel_options, [0 (num-7)/num 1 1/num], obj.owner, 'use_recursive', 'toggle', 'no recursion', 'use recursion');
            obj.input_recursion_depth = GraphicButton(obj.panel_options, [0 (num-8)/num 1 1/num], obj.owner, 'recursion_depth', 'input', 'depth= ');
                        
            obj.button_use_exclude = GraphicButton(obj.panel_options, [0 (num-9)/num 1 1/num], obj.owner, 'use_exclude', 'toggle', 'no exclude', 'use exclude');
            obj.input_exclude_dx = GraphicButton(obj.panel_options, [0 (num-10)/num 1 1/num], obj.owner, 'exclude_dx', 'input', 'exc. dx= ');
            obj.input_exclude_dy = GraphicButton(obj.panel_options, [0 (num-11)/num 1 1/num], obj.owner, 'exclude_dy', 'input', 'exc. dy= ');
                        
            obj.button_use_crop = GraphicButton(obj.panel_options, [0 (num-12)/num 1 1/num], obj.owner, 'use_crop_image', 'toggle', 'no crop', 'use crop');
            obj.input_crop_size = GraphicButton(obj.panel_options, [0 (num-13)/num 1 1/num], obj.owner, 'crop_size', 'input', 'crop= ');
            obj.button_use_save = GraphicButton(obj.panel_options, [0 (num-14)/num 1 1/num], obj.owner, 'use_save_images', 'toggle', 'no save', 'save images');
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%%
            
            obj.panel_info = uipanel('Title', 'info', 'Position', [0.2 0.9 0.8 0.1]);
            
            obj.button_filename = GraphicButton(obj.panel_info, [0.0 0 0.8 1], obj.owner, 'filename', 'info');
%             obj.button_filename.control.HorizontalAlignment = 'right';
            obj.button_filename.font_size = 'small';
            obj.button_filename.Callback = @obj.callback_filename;
            obj.button_batch = GraphicButton(obj.panel_info, [0.8 0 0.1 1], obj.owner, 'batch_num', 'info', 'b= ');
            obj.button_frame = GraphicButton(obj.panel_info, [0.9 0 0.1 1], obj.owner, 'frame_num', 'info', 'f= ');
            
            
            %%%%%%%%%%% panel display %%%%%%%%%%%%%%%%
            
            obj.panel_display = uipanel('Title', 'display', 'Position', [0.0 0.0 1 0.2]);
            
            obj.button_close = GraphicButton(obj.panel_display, [0.0 0.0 0.2 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
            obj.input_line_offset = GraphicButton(obj.panel_display, [0.2 0.5 0.2 0.5], obj.owner, 'line_offset', 'input', ' line offset= ');
            obj.input_rect_size = GraphicButton(obj.panel_display, [0.2 0.0 0.2 0.5], obj.owner, 'rect_size', 'input', ' rect size= ');
            
            obj.button_use_monochrome = GraphicButton(obj.panel_display, [0.4 0.5 0.2 0.5], obj.owner, 'display_monochrome', 'toggle', 'color', 'monochrome');
            
            obj.button_display_which = GraphicButton(obj.panel_display, [0.6 0.5 0.2 0.5], obj.owner, 'cycleDisplayWhich', 'push');
            obj.button_prev = GraphicButton(obj.panel_display, [0.6 0 0.05 0.5], obj.owner, 'prevDisplayIndex', 'push', '-');
            obj.input_streak_num = GraphicButton(obj.panel_display, [0.65 0 0.05 0.5], obj.owner, 'display_index', 'input', ' ');
            obj.button_next = GraphicButton(obj.panel_display, [0.70 0 0.05 0.5], obj.owner, 'nextDisplayIndex', 'push', '+');
            obj.button_best = GraphicButton(obj.panel_display, [0.75 0 0.05 0.5], obj.owner, 'resetDisplayIndex', 'push', 'best');
            
            obj.button_show = GraphicButton(obj.panel_display, [0.8 0.0 0.2 1], obj.owner, 'show', 'push');
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0.2 0.8 0.7]);
            
            obj.update;
            
        end
                            
        function update(obj,~,~)
                        
            if ~obj.check
                return;
            end
           
            for ii = 1:length(obj.buttons)
                obj.buttons{ii}.update;
            end
            
            obj.button_autothresh_val.String = sprintf('%4.1f', obj.owner.latest_thresh);
            
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