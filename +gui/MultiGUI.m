classdef MultiGUI < handle
    
    properties 
        
        owner@radon.MultiFinder; % link back to containg object

        fig@util.plot.FigHandler;
             
        buttons = {};
        
        font_size = 16;
        edit_font_size = 14;
        small_font_size = 12;
        
        debug_bit = 1;
        
    end
    
    properties % gui stuff
        
        panel_info;
        
        panel_controls;
       
        panel_contrast;
        
        panel_objects;
        
        panel_file;
        
        panel_start;
        
        panel_image;
        button_reset_axes;
        axes_image;
        
        panel_close;
        button_close;
        
    end
    
    properties (Hidden=true)
              
        version = 1.00;
        
    end
            
    methods % constructor
       
        function obj = MultiGUI(owner)
            
            % later add other options like copy constructor
            if isa(owner, 'radon.MultiFinder')
                
                if obj.debug_bit, fprintf('MultiGUI constructor v%4.2f\n', obj.version); end
                
                obj.owner = owner;
                
            else
                error('Input a MultiFinder to constructor of MultiGUI!');
            end
            
        end
        
    end
    
    methods % make/update/check GUI
    
        function make(obj)
            
            import util.plot.GraphicButton;
            import util.plot.GraphicPanel;
            import util.plot.ContrastLimits;
            
            obj.buttons = {};
            
            if isempty(obj.fig) || ~isvalid(obj.fig)
                obj.fig = util.plot.FigHandler('MultiFinder');
            end
            
            obj.fig.reset;
            obj.fig.bottom = 5;
            obj.fig.height = 16;
            obj.fig.width = 25;
            obj.fig.name = 'MultiFinder';
            movegui(obj.fig.fig, 'center');
            
            %%%%%%%%%%% panel info %%%%%%%%%%%%%%%%%%
            
            obj.panel_info = GraphicPanel(obj.owner, [0 0.7 0.2 0.3], 'info');
            obj.panel_info.addButton('button_frame', 'frame_index', 'info', 'frame= ', '', '', 0.7);
            obj.panel_info.addButton('button_frame_total', 'total_frames', 'info', '/', '', '', 0.3);
            obj.panel_info.addButton('button_section', 'section_index', 'info', 'section= ', '', '', 0.7);
            obj.panel_info.addButton('button_section_total', 'total_sections', 'info', '/', '', '', 0.3);
            obj.panel_info.addButton('button_thresh', 'current_threshold', 'info', 'thresh= ', '', '', 1);
            obj.panel_info.make;
            
            %%%%%%%%%%% panel controls %%%%%%%%%%%%%%%%%%
            
            obj.panel_controls = GraphicPanel(obj.owner, [0 0.1 0.2 0.6], 'controls');
            obj.panel_controls.number = 9;
%             obj.panel_controls.addButton(
            obj.panel_controls.make;
            
            %%%%%%%%%%% panel contrast %%%%%%%%%%%%%%%%%%
            
            % not sure we need this... 
            
            
            %%%%%%%%%%% panel objects %%%%%%%%%%%%%%%%%%
            
            obj.panel_objects = GraphicPanel(obj.owner, [0.8 0 0.2 1], 'objects');
            obj.panel_objects.number = 10;
            obj.panel_objects.addButton('button_finder', 'finder', 'push', 'Finder');
            obj.panel_objects.addButton('button_classifier', 'class', 'push', 'Classifier');
            obj.panel_objects.make;
        
            %%%%%%%%%%% panel file %%%%%%%%%%%%%%%%%%
            
            obj.panel_file = GraphicPanel(obj.owner, [0.2 0.9 0.6 0.1], 'info');
            obj.panel_file.addButton('button_filename', 'current_filename', 'info', '', '', 'small'); % add option to click this to choose files!
            obj.panel_file.make;
            
            %%%%%%%%%%% panel image %%%%%%%%%%%%%%%%%%
            
            obj.panel_image = uipanel('Title', '', 'Position', [0.2 0.1 0.6 0.8]);
                        
            obj.makeAxes;
            
            obj.button_reset_axes = GraphicButton(obj.panel_image, [0.9 0.95 0.1 0.05], obj.owner, '', 'custom','reset');
            obj.button_reset_axes.Callback = @obj.makeAxes;
            
            %%%%%%%%%%% panel start %%%%%%%%%%%%%%%%%%
            
            obj.panel_start = GraphicPanel(obj.owner, [0.2 0 0.8 0.1], '');
            obj.panel_start.addButton('button_start', '', 'custom', 'RUN', '', '', 0.8);
            obj.panel_start.addButton('button_reset', 'reset', 'push', 'RESET', '', '', 0.2);
            obj.panel_start.make;
            obj.panel_start.button_reset.Callback = @obj.callback_run;
            
            %%%%%%%%%%% panel close %%%%%%%%%%%%%%%%%%
            
            obj.panel_close = uipanel('Position', [0 0 0.2 0.1]);
            obj.button_close = GraphicButton(obj.panel_close, [0 0 1 1], obj.owner, '', 'custom', 'CLOSE');
            obj.button_close.Callback = @obj.callback_close;
            
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
            
            if obj.owner.brake_bit
                obj.panel_start.button_start.String = 'RUN';
            else
                obj.panel_start.button_start.String = 'STOP';
            end
            
            
        end
                        
        function c = check(obj)
           
            c = ~isempty(obj) && ~isempty(obj.panel_close) && isvalid(obj.panel_close);
            
        end
        
    end
                
    methods % callbacks
        
        function callback_run(obj, ~, ~)
           
            
            if obj.owner.brake_bit
                if obj.debug_bit, disp('callback: run'); end            
                obj.owner.run;
            else
                if obj.debug_bit, disp('callback: stop'); end            
                obj.owner.brake_bit = 1;
            end
            
            obj.update;
            
        end
        
        function callback_close(obj, ~, ~)
           
            if obj.debug_bit, disp('callback: close'); end
            
            delete(obj.fig.fig);
            
        end
        
    end
    
end