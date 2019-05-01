classdef ImageComponentParser < hgsetget
    % class ImageComponentParser: gui to analyze the component
    % decompositions of imaging data.
    %
    % @file: ImageComponentParser.m
    % @author: Paxon Frady
    % @created: 8/20/2013
    
    properties
        % gui properties
        h; % graphic object handles.
        gui; % settings for the gui.
        
        % hgsetget properties
        Position; % The size and position of the object
        Parent; % The parent of the object
        
        % object properties
        im_data; % The original image data

        rois; % Cell array containing the sets of rois
        
        settings; % struct containing settings for each stage.
        pp;  % struct containing preprocessing data
        pca; % struct containing pca data
        ica; % struct containing ica data
        post; % struct containing postprocessing data
        cluster; % struct containing clustering data
        
        stage; % Keeps track of which stage the analysis is in.
    end
    
    properties (Constant)
        % The stage values need to go in increasing order
        InitStage = 0;
        PreProcessingStage = 1;
        PcaStage = 2;
        IcaStage = 3;
        PostProcessingStage = 4;
    end
    
    events
        
    end
    
    methods
        %%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = ImageComponentParser(parent, im_data)
            % self = ImageComponentParser(parent, im_data): parses the
            % components of imaging data.
            %
            % @param: parent the parent handle of the gui object. If no
            % parent is given a new figure will be created and used as the
            % parent.
            % @param: im_data MxNxt matrix of image data.
            % @return: self handle to the ImageComponentParser instance
            
            if nargin < 1
                parent = [];
            end
            if nargin < 2
                im_data = rand(128, 128, 10);
            end
            
            self.im_data = im_data;
            
            self = self.init_state();
            self = self.init_gui(parent);
            
            self.update();
        end
        
        function self = init_state(self)
            % init_state: initializes the gui state variables.
            
            self.Position = [0 0 1024 700];
            
            self.default_settings();
            
            self.rois{1} = [];
            self.gui.current_roi_set = 1;
            self.gui.roi_names{1} = 'ROI_set_1';
            self.gui.is_component_set = false;

            
            self.gui.MARGIN = 10;
            self.gui.CAX_H = 240;
            self.gui.BUTTON_H = 30;
            self.gui.PC_PANEL_ITEM_W = 20;
            
            % 3d-axis state variables
            self.gui.x_pc = 1; % The PC scores to plot on the x-axis
            self.gui.y_pc = 2; % The PC scores to plot on the y-axis
            self.gui.z_pc = 3; % The PC scores to plot on the z-axis 
            
            self.gui.rotate_rate = 0.05;
            self.gui.rotate_angle_step = 2;
            self.gui.rotate_timer = timer('TimerFcn', @self.rotate_timer_cb, ...
                                          'ExecutionMode', 'fixedRate', ...
                                          'Period', self.gui.rotate_rate);
                                      
            self.gui.locked_components = {};
            self.gui.locked_colors = lines(7);
            
            self.gui.display = 1;
            
            self.gui.data_norm_frame = 20; % Data normalized to this frame
            
            self.stage = self.InitStage;
        end
        
        function self = init_gui(self, parent)
            % init_gui: initializes the gui objects.
            %
            % @param: parent the parent handle of this object. Default is
            % to create a new figure. Use [] for default.
            
            if nargin < 2 || isempty(parent)
                % No parent given, then make a new figure as the parent.
                parent = figure();
                clf(parent);
                set(parent, 'Position', [100, 100, self.Position(3), self.Position(4)]);
            end
            
            self.h.Parent = parent;
            
            % We must use a figure event notifier, which is attached to the
            % top-most figure. Go up until we get the figure.
            self.h.fh = self.h.Parent;
            while ~strcmp(get(self.h.fh, 'Type'), 'figure')
                % Then the fh object is not a figure, so go up.
                self.h.fh = get(self.h.fh, 'Parent');
            end
            
            % Now fh is the parent figure. Set its event notifier.
            self.gui.fen = FigureEventNotifier(self.h.fh);
            addlistener(self.gui.fen, 'WindowKeyPress', @self.window_key_press_cb);
            set(self.h.fh, 'Toolbar', 'figure');
            
            %self.h.panel = uipanel(self.Parent, 'Units', 'pixels');
            self.h.panel = uiextras.BoxPanel('Parent', self.h.Parent, ...
                'Title', 'Image Component Parser');
            
            % RoiEditor
            self.h.roi_editor = RoiEditor([], self.im_data);
            %self.h.roi_editor.is_editing_enabled = false;
            
            addlistener(self.h.roi_editor, 'SelectionChanged', @self.re_selection_changed_cb);
            addlistener(self.h.roi_editor, 'NewRois', @self.re_new_rois_cb);
            addlistener(self.h.roi_editor, 'AlteredRois', @self.re_altered_rois_cb);
            addlistener(self.h.roi_editor, 'DeletedRois', @self.re_deleted_rois_cb);
            %%%
            addlistener(self.h.roi_editor, 'FrameChanged', @self.re_frame_changed_cb);
            
            % Roi Management
            self.h.roi_listbox = uicontrol('Style', 'listbox', ...
                'String', self.gui.roi_names, 'Min', 1, 'Max', 1, 'Callback', @self.roi_listbox_cb);
            self.h.add_rois_button = uicontrol('Style', 'pushbutton', ...
                'String', '+', 'Callback', @self.add_rois_button_cb);
            self.h.delete_rois_button = uicontrol('Style', 'pushbutton', ...
                'String', '-', 'Callback', @self.delete_rois_button_cb);
            
            % Component axes
            self.h.component_axes = axes('Units', 'pixels');
            self.h.lock_component_button = uicontrol('Style', 'pushbutton',...
                'String', 'Lock', 'Callback', @self.lock_component_button_cb);
            self.h.clear_locked_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Clear', 'Callback', @self.clear_locked_button_cb);
            
            % Data axes
            self.h.data_axes = axes('Units', 'pixels');
                        
             
            self.h.stage_tab_panel = uiextras.TabPanel('Callback', @self.stage_tab_panel_cb);
            self.h.data_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.pp_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.pca_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.ica_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.post_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.cluster_tab = uiextras.Panel('Parent', self.h.stage_tab_panel);
            self.h.stage_tab_panel.TabNames = {'Data', 'Pre-Processing', 'PCA', 'ICA', 'Segment', 'Cluster'};
            self.h.stage_tab_panel.SelectedChild = self.gui.display;
            
            self.h.stage_status = uicontrol('Style', 'text', 'String', '...');
            
            % Data tab
            self.h.reset_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Reset', 'Callback', @self.reset_button_cb);
            
            % PP tab
            self.h.motion_correct_button = uicontrol('Style', 'pushbutton',...
                'String', 'Motion Correct', 'Callback', @self.motion_correct_button_cb);
            self.h.preprocess_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Pre Process', 'Callback', @self.preprocess_button_cb);
            self.h.wavelets_check = uicontrol('Style', 'checkbox', 'Value', 0, ...
                'String', 'Use Wavelets', 'Callback', @self.wavelets_check_cb);
            
            % PCA tab
            self.h.runpca_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Run PCA', 'Callback', @self.runpca_button_cb);            
            
            % ICA tab
            self.h.runica_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Run ICA', 'Callback', @self.runica_button_cb);

            % Post tab
            self.h.calc_rois_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Calc ROIs', 'Callback', @self.calc_rois_cb);
            self.h.segment_ics_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Segment ICs', 'Callback', @self.segment_ics_cb);
            
            % Cluster tab
            self.h.estimate_clusters_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Estimate Clusters', 'Callback', @self.estimate_clusters_cb);
            self.h.set_cluster_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Set Cluster', 'Callback', @self.set_cluster_button_cb);
            self.h.uncluster_button = uicontrol('Style', 'pushbutton', ...
                'String', 'Uncluster', 'Callback', @self.uncluster_button_cb);
            
            
            % Axes to display three PCs
            self.h.pc3_axes = axes('Units', 'pixels');  
            % Popups for the plot.
            self.h.x_popup = uicontrol('Style', 'popupmenu', ...
                'Units', 'pixels', 'String', {'None'});
            self.h.y_popup = uicontrol('Style', 'popupmenu', ...
                'Units', 'pixels', 'String', {'None'});
            self.h.z_popup = uicontrol('Style', 'popupmenu', ...
                'Units', 'pixels', 'String', {'None'});
            set(self.h.x_popup, 'Callback', @self.x_popup_cb);
            set(self.h.y_popup, 'Callback', @self.y_popup_cb);
            set(self.h.z_popup, 'Callback', @self.z_popup_cb);
            % Buttons to auto-rotate to the 2-d plots.
            self.h.xy_button = uicontrol('Style', 'pushbutton', ...
                'Units', 'pixels', 'String', 'xy', 'Callback', @self.xy_button_cb);
            self.h.yz_button = uicontrol('Style', 'pushbutton', ...
                'Units', 'pixels', 'String', 'yz', 'Callback', @self.yz_button_cb);
            self.h.zx_button = uicontrol('Style', 'pushbutton', ...
                'Units', 'pixels', 'String', 'xz', 'Callback', @self.zx_button_cb);
            % Labels for the plot controls. Even though these are labels,
            % strange things were happening when I didn't keep the handles.
            self.h.xl = uicontrol('Style', 'text', 'String', 'X:');
            self.h.yl = uicontrol('Style', 'text', 'String', 'Y:');
            self.h.zl = uicontrol('Style', 'text', 'String', 'Z:');
            % Toggle button to turn continuous rotation on and off.
            self.h.rotate_toggle = uicontrol('Style', 'togglebutton', ...
                'Units', 'pixels', 'String', 'Rotate', 'Callback', @self.rotate_toggle_cb);
            
            % Menu
            self.h.main_menu = uimenu(self.h.fh, 'Label', 'ICP');
            self.h.load_settings_item = uimenu('Parent', self.h.main_menu, ...
                'Label', 'Load Settings', 'Callback', @self.load_settings_cb);
            self.h.save_settings_item = uimenu('Parent', self.h.main_menu, ...
                'Label', 'Save Settings', 'Callback', @self.save_settings_cb);
            set(self.h.save_settings_item, 'Separator', 'on');
            self.h.edit_settings_item = uimenu('Parent', self.h.main_menu, ...
                'Label', 'Edit Settings');
            self.h.edit_preprocessing_item = uimenu('Parent', self.h.edit_settings_item, ...
                'Label', 'Preprocessing', 'Callback', @self.edit_preprocessing_cb);
            self.h.edit_pca_item = uimenu('Parent', self.h.edit_settings_item, ...
                'Label', 'PCA', 'Callback', @self.edit_pca_cb);            
            self.h.edit_ica_item = uimenu('Parent', self.h.edit_settings_item, ...
                'Label', 'ICA', 'Callback', @self.edit_ica_cb);
            
            self.h.save_rois_item = uimenu('Parent', self.h.main_menu, ...
                'Label', 'Save ROIs', 'Callback', @self.save_rois_cb,...
                'Separator', 'on');
            self.h.load_rois_item = uimenu('Parent', self.h.main_menu, ...
                'Label', 'Load ROIs', 'Callback', @self.load_rois_cb);
            
            %self = self.reset_layout();
            self = self.uiextras_layout();
        end
        
        function self = uiextras_layout(self)
            % reset_layout: Sets the layout of all of the gui components.
            
            % Ok, I'm going to do this with the layouts, and reset the
            % parents here.
            disp('uiextras_layout');
            % 3D plot Layout elements
            self.h.pc3_main_panel = uiextras.BoxPanel('Title', '3D Plot');
            self.h.pc3_main_hbox = uiextras.HBox();
            self.h.pc3_control_rotate_vbox = uiextras.VBox();
            self.h.pc3_control_panel = uiextras.BoxPanel('Title', '3D Plot Controls');
            self.h.pc3_control_grid = uiextras.Grid();
            
            % Component axes elements
            self.h.component_panel = uiextras.BoxPanel('Title', 'Component Plot');
            self.h.component_lock_vbox = uiextras.VBox();
            self.h.lock_button_hbox = uiextras.HButtonBox();
            self.h.component_pc3_vbox = uiextras.VBoxFlex();
            
            % Data axes elements
            self.h.data_panel = uiextras.BoxPanel('Title', 'Data Plot');
            
            % Stage elements
            %self.h.stage_button_hbox = uiextras.HButtonBox();
            self.h.editor_button_vbox = uiextras.VBox();   
            
            self.h.roi_stage_hbox = uiextras.HBox();
            self.h.roi_control_vbox = uiextras.VBox();
            self.h.roi_button_hbox = uiextras.HButtonBox();
            
            self.h.cluster_button_hbox = uiextras.HButtonBox();
            self.h.post_button_hbox = uiextras.HButtonBox();
            self.h.ica_button_hbox = uiextras.HButtonBox();
            self.h.pca_button_hbox = uiextras.HButtonBox();
            self.h.pp_button_hbox = uiextras.HButtonBox();
            self.h.data_button_hbox = uiextras.HButtonBox();
            
            self.h.main_hbox = uiextras.HBoxFlex();
            
            % Top-level hierarchy
            set(self.h.main_hbox, 'Parent', self.h.panel);
            set(self.h.editor_button_vbox, 'Parent', self.h.main_hbox);
            set(self.h.component_pc3_vbox, 'Parent', self.h.main_hbox);

            % Button, ROI Editor hierarchy
            set(self.h.roi_editor, 'Parent', self.h.editor_button_vbox);
            set(self.h.roi_stage_hbox, 'Parent', self.h.editor_button_vbox);
            set(self.h.stage_status, 'Parent', self.h.editor_button_vbox.double());
            
            set(self.h.roi_control_vbox, 'Parent', self.h.roi_stage_hbox);
            set(self.h.stage_tab_panel, 'Parent', self.h.roi_stage_hbox);
            
            set(self.h.roi_listbox, 'Parent', self.h.roi_control_vbox.double());
            set(self.h.roi_button_hbox, 'Parent', self.h.roi_control_vbox);
            
            set(self.h.add_rois_button, 'Parent', self.h.roi_button_hbox.double());
            set(self.h.delete_rois_button, 'Parent', self.h.roi_button_hbox.double());            
            
            % Data tab
            set(self.h.data_button_hbox, 'Parent', self.h.data_tab);
            set(self.h.reset_button, 'Parent', self.h.data_button_hbox.double());
            
            % PP tab
            set(self.h.pp_button_hbox, 'Parent', self.h.pp_tab);
            set(self.h.motion_correct_button, 'Parent', self.h.pp_button_hbox.double());
            set(self.h.wavelets_check, 'Parent', self.h.pp_button_hbox.double());
            set(self.h.preprocess_button, 'Parent', self.h.pp_button_hbox.double());            
            
            % PCA tab
            set(self.h.pca_button_hbox, 'Parent', self.h.pca_tab);
            set(self.h.runpca_button, 'Parent', self.h.pca_button_hbox.double());
            
            % ICA tab
            set(self.h.ica_button_hbox, 'Parent', self.h.ica_tab)
            set(self.h.runica_button, 'Parent', self.h.ica_button_hbox.double());
            
            % Post tab
            set(self.h.post_button_hbox, 'Parent', self.h.post_tab);
            set(self.h.calc_rois_button, 'Parent', self.h.post_button_hbox.double());
            set(self.h.segment_ics_button, 'Parent', self.h.post_button_hbox.double());
            
            % Cluster tab
            set(self.h.cluster_button_hbox, 'Parent', self.h.cluster_tab);
            set(self.h.estimate_clusters_button, 'Parent', self.h.cluster_button_hbox.double());
            set(self.h.set_cluster_button, 'Parent', self.h.cluster_button_hbox.double());
            set(self.h.uncluster_button, 'Parent', self.h.cluster_button_hbox.double());
            
            % Component axis hierarchy
            set(self.h.data_panel, 'Parent', self.h.component_pc3_vbox);
            set(self.h.component_panel, 'Parent', self.h.component_pc3_vbox);
            set(self.h.pc3_main_panel, 'Parent', self.h.component_pc3_vbox);
            
            set(self.h.component_lock_vbox, 'Parent', self.h.component_panel);
            set(self.h.component_axes, 'Parent', self.h.component_lock_vbox.double());
            set(self.h.lock_button_hbox, 'Parent', self.h.component_lock_vbox);
            set(self.h.lock_component_button, 'Parent', self.h.lock_button_hbox.double());
            set(self.h.clear_locked_button, 'Parent', self.h.lock_button_hbox.double());
            
            set(self.h.data_axes, 'Parent', self.h.data_panel.double());
            
            % 3D plot hierarchy
            set(self.h.pc3_main_hbox, 'Parent', self.h.pc3_main_panel);
            
            set(self.h.pc3_axes, 'Parent', self.h.pc3_main_hbox.double());
            set(self.h.pc3_control_rotate_vbox, 'Parent', self.h.pc3_main_hbox);
            
            set(self.h.pc3_control_panel, 'Parent', self.h.pc3_control_rotate_vbox);
            set(self.h.rotate_toggle, 'Parent', self.h.pc3_control_rotate_vbox.double());
            
            set(self.h.pc3_control_grid, 'Parent', self.h.pc3_control_panel);
                 
            set(self.h.xl, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.yl, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.zl, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.x_popup, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.y_popup, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.z_popup, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.xy_button, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.yz_button, 'Parent', self.h.pc3_control_grid.double());
            set(self.h.zx_button, 'Parent', self.h.pc3_control_grid.double());
            
            % Layout structure.
            set(self.h.pc3_control_grid,...
                'ColumnSizes', [self.gui.PC_PANEL_ITEM_W, -1, self.gui.PC_PANEL_ITEM_W],...
                'RowSizes', [self.gui.BUTTON_H, self.gui.BUTTON_H, self.gui.BUTTON_H]);
            set(self.h.pc3_control_rotate_vbox, 'Sizes', [-1, self.gui.BUTTON_H]);
            set(self.h.pc3_main_hbox, 'Sizes', [-4, -1]);
            
            set(self.h.component_pc3_vbox, 'Sizes', [-1, -1, -1]);
            
            set(self.h.component_lock_vbox, 'Sizes', [-1, self.gui.BUTTON_H]);
            
            set(self.h.editor_button_vbox, 'Sizes', [550, -1, self.gui.BUTTON_H]);
            set(self.h.roi_stage_hbox, 'Sizes', [-1, -3]);
            set(self.h.roi_control_vbox, 'Sizes', [-1, self.gui.BUTTON_H]);
            
            set(self.h.main_hbox, 'Sizes', [515 -1]);
            
            set(self.h.Parent, 'Position', [100, 100, self.Position(3), self.Position(4)]);
        end
        
        function self = reset_layout(self)
            % reset_layout: Sets the layout of all of the gui components.
            
%             set(self.h.panel, 'Position', self.Position);
% 
%             PANEL_W = self.Position(3);
%             PANEL_H = self.Position(4);      
%             
%             %%% ROI Editor
%             RE_L = self.gui.MARGIN;
%             RE_B = self.gui.MARGIN + 2*self.gui.BUTTON_H; %PANEL_H - RE_H - self.gui.MARGIN;
%             RE_W = self.h.roi_editor.Position(3);
%             RE_H = PANEL_H - RE_B - self.gui.MARGIN; %self.h.roi_editor.Position(4);
%             
%             self.h.roi_editor.Position = [RE_L, RE_B, RE_W, RE_H];
%             self.h.roi_editor.reset_layout();
%             
%             %%% Stage Buttons
%             B_W = RE_W / 4;
%             B_H = self.gui.BUTTON_H;
%             B_L = @(N) N * B_W + self.gui.MARGIN;
%             B_B = self.gui.MARGIN + self.gui.BUTTON_H;
%             
%             set(self.h.reset_button, 'Position', [B_L(0), B_B, B_W, B_H]);
%             set(self.h.preprocess_button, 'Position', [B_L(1), B_B, B_W, B_H]);
%             set(self.h.runpca_button, 'Position', [B_L(2), B_B, B_W, B_H]);
%             set(self.h.runica_button, 'Position', [B_L(3), B_B, B_W, B_H]);
%             
%             %%% Component Axes
%             CAX_L = 3 * self.gui.MARGIN + RE_W;
%             CAX_W = PANEL_W - CAX_L - self.gui.MARGIN;
%             CAX_H = self.gui.CAX_H;
%             CAX_B = PANEL_H - CAX_H - self.gui.MARGIN;
%             set(self.h.component_axes, 'Position', [CAX_L, CAX_B, CAX_W, CAX_H]);
% 
%             %%% Score Axes
%             SAX_L = CAX_L;
%             SAX_B = RE_B;
%             SAX_W = CAX_W - self.gui.PC_PANEL_W;
%             SAX_H = CAX_B - SAX_B - 3 * self.gui.MARGIN;
%             set(self.h.pc3_axes, 'Position', [SAX_L, SAX_B, SAX_W, SAX_H]);
%             axis(self.h.pc3_axes, 'vis3d');
%             
%             %%% Score Controls
%             LABEL_W = 20;
%             PUP_W = self.gui.PC_PANEL_W - 2 * LABEL_W - self.gui.MARGIN;
%             PUP_L = PANEL_W - self.gui.PC_PANEL_W;
%             PUP_B = @(N) SAX_B + SAX_H - self.gui.MARGIN - (N + 1) * self.gui.BUTTON_H;
%             PUP_H = self.gui.BUTTON_H;
%             
%             set(self.h.xl, 'Position', [PUP_L, PUP_B(0), LABEL_W, PUP_H]);
%             set(self.h.x_popup, 'Position', [PUP_L + LABEL_W, PUP_B(0), PUP_W, PUP_H]);
%             set(self.h.xy_button, 'Position', [PUP_L + LABEL_W + PUP_W, PUP_B(0), LABEL_W, PUP_H]);
%             
%             set(self.h.yl, 'Position', [PUP_L, PUP_B(1), LABEL_W, PUP_H]);
%             set(self.h.y_popup, 'Position', [PUP_L + LABEL_W, PUP_B(1), PUP_W, PUP_H]);
%             set(self.h.yz_button, 'Position', [PUP_L + LABEL_W + PUP_W, PUP_B(1), LABEL_W, PUP_H]);
%             
%             set(self.h.zl, 'Position', [PUP_L, PUP_B(2), LABEL_W, PUP_H]);
%             set(self.h.z_popup, 'Position', [PUP_L + LABEL_W, PUP_B(2), PUP_W, PUP_H]);
%             set(self.h.zx_button, 'Position', [PUP_L + LABEL_W + PUP_W, PUP_B(2), LABEL_W, PUP_H]);
%             
%             set(self.h.rotate_toggle, 'Position', [PUP_L, PUP_B(4), self.gui.PC_PANEL_W, PUP_H]);
        end
        
        %%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function window_key_press_cb(self, source_h, eventdata)
            disp('window_keypress_cb');
        end
        
        function reset_button_cb(self, source_h, eventdata)
            disp('reset_button_cb');
            self.stage = self.InitStage;
            self.update();
        end
        
        function motion_correct_button_cb(self, source_h, eventdata)
            disp('motion_correct_button_cb');
            self.run_motion_correction();
        end
        
        function preprocess_button_cb(self, source_h, eventdata)
            disp('preprocess_button_cb');
            
            self.run_preprocessing();
        end
        
        function wavelets_check_cb(self, source_h, eventdata)
            disp('wavelets_check_cb');
            disp(get(source_h, 'Value'));
        end
        
        function runpca_button_cb(self, source_h, eventdata)
            disp('runpca_button_cb');
            
            self.run_pca();
        end
        
        function runica_button_cb(self, source_h, eventdata)
            disp('runica_button_cb');
            
            self.run_ica();
        end
        
        function segment_ics_cb(self, source_h, eventdata)
            disp('segment_ics_cb');
            
            self.run_segmentation();
        end
        
        function estimate_clusters_cb(self, source_h, eventdata)
            self.estimate_clusters();
        end
        
        function set_cluster_button_cb(self, source_h, eventdata)
            disp('set_cluster_button_cb');
            
            self.set_cluster();
        end
        
        function uncluster_button_cb(self, source_h, eventdata)
            disp('uncluster_button_cb');
            
            self.uncluster();
        end
        
        function stage_tab_panel_cb(self, source_h, eventdata)
            disp('stage_tab_panel_cb');
            
            self.gui.display = eventdata.SelectedChild;
            self.update();
        end
        
        function x_popup_cb(self, sh, ed)
            % x_popup_cb: callback for when x popup is changed.
            
            v = get(sh, 'Value');
            self.gui.x_pc = v;
            
            % Now make sure y, z aren't x.
            if self.gui.y_pc == self.gui.x_pc
                % Then y is equal to x
                for i = 1:10
                    if i ~= self.gui.x_pc && i ~= self.gui.z_pc
                        self.gui.y_pc = i;
                        break;
                    end 
                end
            end
            if self.gui.z_pc == self.gui.x_pc
                % Then z is equal to x
                for i = 1:10
                    if i ~= self.gui.x_pc && i ~= self.gui.y_pc
                        self.gui.z_pc = i;
                        break;
                    end 
                end
            end
            self.update();
        end
        
        function y_popup_cb(self, sh, ed)
            % y_popup_cb: callback for when y popup is changed.
            
            v = get(sh, 'Value');
            self.gui.y_pc = v;
            
            % Now make sure x, z aren't y.
            if self.gui.x_pc == self.gui.y_pc
                % Then x is equal to y
                for i = 1:10
                    if i ~= self.gui.y_pc && i ~= self.gui.z_pc
                        self.gui.x_pc = i;
                        break;
                    end 
                end
            end
            if self.gui.z_pc == self.gui.y_pc
                % Then z is equal to y
                for i = 1:10
                    if i ~= self.gui.y_pc && i ~= self.gui.x_pc
                        self.gui.z_pc = i;
                        break;
                    end 
                end
            end
            self.update();
        end
        
        function z_popup_cb(self, sh, ed)
            % z_popup_cb: callback for when z popup is changed.
            
            v = get(sh, 'Value');
            self.gui.z_pc = v;
            
            % Now make sure y, x aren't z.
            if self.gui.y_pc == self.gui.z_pc
                % Then y is equal to z
                for i = 1:10
                    if i ~= self.gui.x_pc && i ~= self.gui.z_pc
                        self.gui.y_pc = i;
                        break;
                    end 
                end
            end
            if self.gui.x_pc == self.gui.z_pc
                % Then x is equal to z
                for i = 1:10
                    if i ~= self.gui.z_pc && i ~= self.gui.y_pc
                        self.gui.x_pc = i;
                        break;
                    end 
                end
            end
            self.update();
        end
        
        function xy_button_cb(self, sh, ed)
            % xy_button_cb: callback that rotates pc axes to show xy.
            view(self.h.pc3_axes, [0 90]);
        end
        
        function yz_button_cb(self, sh, ed)
            % yz_button_cb: callback that rotates pc axes to show yz.
            view(self.h.pc3_axes, [90 0]);
        end
        
        function zx_button_cb(self, sh, ed)
            % zx_button_cb: callback that rotates pc axes to show zx.
            view(self.h.pc3_axes, [0 0]);
        end
        
        function rotate_toggle_cb(self, sh, ed)
            v = get(sh, 'Value');
            
            if v
                start(self.gui.rotate_timer);
            else
                stop(self.gui.rotate_timer);
            end
        end
        
        function rotate_timer_cb(self, sh, ed)
            [az, el] = view(self.h.pc3_axes);
            view(self.h.pc3_axes, az + self.gui.rotate_angle_step, el);
        end
        
        function re_selection_changed_cb(self, source_h, eventdata)
            disp('re_selection_changed_cb');
            self.update();
        end
        
        function re_new_rois_cb(self, source_h, eventdata)
            disp('re_new_rois_cb');
            self.update_current_roi_set();
            self.update();
        end
        
        function re_altered_rois_cb(self, source_h, eventdata)
            disp('re_altered_rois_cb');
            self.update_current_roi_set();
            self.update();
        end
        
        function re_deleted_rois_cb(self, source_h, eventdata)
            disp('re_deleted_rois_cb');
            self.update_current_roi_set();
            self.update();
        end
        
        function re_frame_changed_cb(self, source_h, eventdata)
            disp('re_frame_changed_cb');
            self.update();
        end
        
        function load_settings_cb(self, source_h, eventdata)
            disp('load_settings_cb');
        end
        
        function save_settings_cb(self, source_h, eventdata)
            disp('save_settings_cb');
        end
        
        function edit_preprocessing_cb(self, source_h, eventdata)
            disp('edit_preprocessing_cb');
            
            self.edit_preprocessing_settings();
        end
        
        function edit_pca_cb(self, source_h, eventdata)
            disp('edit_pca_cb');
        end
        
        function edit_ica_cb(self, source_h, eventdata)
            disp('edit_ica_cb');
        end
        
        function pp_settings_ok_cb(self, source_h, eventdata)
            
            % Get all of the data from the dialog. If the dialog has bad
            % values, then these will be NaN. The set function will handle
            % the NaN case.
            smooth_M_val = str2double(get(self.h.pp_smooth_M, 'String'));
            smooth_N_val = str2double(get(self.h.pp_smooth_N, 'String'));
            smooth_T_val = str2double(get(self.h.pp_smooth_T, 'String'));
            
            down_M_val = str2double(get(self.h.pp_down_M, 'String'));
            down_N_val = str2double(get(self.h.pp_down_N, 'String'));
            down_T_val = str2double(get(self.h.pp_down_T, 'String'));
            
            self.set_pp_smooth_window(smooth_M_val, smooth_N_val, smooth_T_val);
            self.set_pp_down_sample(down_M_val, down_N_val, down_T_val);
            
            close(self.h.pp_dialog);
        end
        
        function pp_settings_cancel_cb(self, source_h, eventdata)
            close(self.h.pp_dialog);
        end
        
        function lock_component_button_cb(self, source_h, eventdata)
            disp('lock_component_button_cb');
            self.lock_current_component();
        end
        
        function clear_locked_button_cb(self, source_h, eventdata)
            disp('clear_locked_button_cb');
            self.clear_locked_components();
        end
        
        function component_label_cb(self, source_h, eventdata)
            disp('component_label_cb');
            c = get(source_h, 'UserData');
            
            if isscalar(c)
                self.h.roi_editor.set_current_frame(c);
            end
        end
        
        function roi_listbox_cb(self, source_h, eventdata)
            disp('roi_listbox_cb');
            
            sel = get(self.h.roi_listbox, 'Value');
            self.select_roi_set(sel);
        end
        
        function add_rois_button_cb(self, source_h, eventdata)
            disp('add_rois_button_cb');
            self.create_new_roi_set();
        end
        
        function delete_rois_button_cb(self, source_h, eventdata)
            disp('delete_rois_button_cb');
            self.delete_roi_set();
        end
        
        function calc_rois_cb(self, source_h, eventdata)
            disp('calc_rois_cb');
            self.calc_rois();
        end
        
        function save_rois_cb(self, source_h, eventdata)
            self.save_rois();
        end
        
        function load_rois_cb(self, source_h, eventdata)
            self.load_rois();
        end
        
        %%%%%%%%%% Main Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function update(self)
            % update(self): main update function
            
            disp('ICP: update');
            f = self.h.roi_editor.current_frame;
            
            rois = self.h.roi_editor.rois.xyrra(self.h.roi_editor.selected_rois, :, 1);
            
            cla(self.h.data_axes);
            if ~isempty(rois)
                % plot the original data from the ROI.
                data = mean_roi(self.im_data, rois);
                
                data = 100 * (data ./ repmat(data(self.gui.data_norm_frame, :), size(data, 1), 1) - 1);
                
                plot(self.h.data_axes, data);
            end
            
            if self.gui.display == 6 && self.stage >= self.PostProcessingStage ...
                    && isfield(self.cluster, 'im_data') && ~isempty(self.cluster.im_data)
                self.h.roi_editor.set_im_data(self.cluster.im_data, [], self.cluster.im_x, self.cluster.im_y);
                f = self.h.roi_editor.set_current_frame(f);
            
                cla(self.h.component_axes);
                plot(self.h.component_axes, self.ica.ics(:, f), 'k', 'LineWidth', 2);
                
                plot(self.h.component_axes, self.ica.ics(:, self.cluster.top3(f, 3)), 'b');
                plot(self.h.component_axes, self.ica.ics(:, self.cluster.top3(f, 2)), 'g');
                plot(self.h.component_axes, self.ica.ics(:, self.cluster.top3(f, 1)), 'r');
                
            elseif self.gui.display == 5 && self.stage >= self.PostProcessingStage ...
                    && isfield(self.post, 'im_data') && ~isempty(self.post.im_data)
                % Then view the segments
                self.h.roi_editor.set_im_data(self.post.im_data, [], self.post.im_x, self.post.im_y);
                f = self.h.roi_editor.set_current_frame(f);
                
            elseif self.gui.display == 4 && self.stage >= self.IcaStage
                % Then view the ICs
                self.h.roi_editor.set_im_data(self.ica.im_data, [], self.ica.im_x, self.ica.im_y);
                f = self.h.roi_editor.set_current_frame(f);
                
                % Plot the current ic
                cla(self.h.component_axes);
                legend(self.h.component_axes, 'off');
                
                plot(self.h.component_axes, self.ica.ics(:, f), 'k');
                
            elseif self.gui.display == 3 && self.stage >= self.PcaStage
                % Then view the PCs
                self.h.roi_editor.set_im_data(self.pca.im_data, [], self.pca.im_x, self.pca.im_y);
                f = self.h.roi_editor.set_current_frame(f);
                
                % Plot the current pc
                cla(self.h.component_axes);
                legend(self.h.component_axes, 'off');
                
                plot(self.h.component_axes, self.pca.pcs(:, f), 'k');
            elseif self.gui.display == 2 && self.stage >= self.PreProcessingStage
                % View the Preprocessing data
                self.h.roi_editor.set_im_data(self.pp.im_data, [], self.pp.im_x, self.pp.im_y);
                self.h.roi_editor.set_current_frame(f);
            else
                % View the raw data. 
                self.h.roi_editor.set_im_data(self.im_data);
                self.h.roi_editor.set_current_frame(f);
            end
            
            switch self.stage
                case self.InitStage, s = 'Initialized';
                case self.PreProcessingStage, s = 'Pre-Processing Complete';
                case self.PcaStage, s = 'PCA Complete';
                case self.IcaStage, s = 'ICA Complete';
                case self.PostProcessingStage, s = 'Post-Processing Complete';
            end
                    
            set(self.h.stage_status, 'String', s);
            
            
            self.h.roi_editor.color_rois = true;
            if self.gui.display == 6 && self.stage == self.PostProcessingStage
                self.h.roi_editor.color_rois = false;
                self.draw_clustered_rois();
                
            elseif self.gui.is_component_set(self.gui.current_roi_set) ...
                    && self.stage >= self.IcaStage ...
                    && ~isempty(self.h.roi_editor.selected_rois)
                cla(self.h.component_axes);
                
                self.h.component_labels = [];
                
                for i = 1:length(self.h.roi_editor.selected_rois)
                    c = self.gui.locked_colors(mod(i-1, length(self.gui.locked_colors))+1, :);
                    
                    plot(self.h.component_axes, self.ica.ics(:, self.h.roi_editor.selected_rois(i)), 'Color', c);
                    
                    text_x = 0.95 - (length(self.h.roi_editor.selected_rois) - i) * 0.05;
                    self.h.component_labels(i) = text(text_x, 0.9, num2str(self.h.roi_editor.selected_rois(i)), ...
                        'Parent', self.h.component_axes, 'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', 'UserData', self.h.roi_editor.selected_rois(i), ...
                        'Color', c, 'ButtonDownFcn', @self.component_label_cb);
                end
            else
                self.show_best_selected_components(rois);
            end
            
            self.plot_locked_components();
            self.plot_3d();
            
            self.build_roi_listbox();
            
            drawnow;
        end
        
        function run_motion_correction(self)
            % run_motion_correction(self): Runs the motion correction.

            disp('run_motion_correction');
            
            set(self.h.stage_status, 'String', 'Running Motion Correction...');
            drawnow;
            
            s = self.settings.preprocessing;
            
            [self.pp.im_data_mc, warp] = s.motion_correct_func(self.im_data);
            
            % Fix the size of the corrected data.
            sx = floor((size(self.im_data, 2) - size(self.pp.im_data_mc, 2))/2) + 1; 
            sy = floor((size(self.im_data, 1) - size(self.pp.im_data_mc, 1))/2) + 1;
            
            self.pp.mc_x = sx:(size(self.pp.im_data_mc, 2) + sx - 1);
            self.pp.mc_y = sy:(size(self.pp.im_data_mc, 1) + sy - 1);
        end
        
        function run_filters(self)
            % run_filters(self): Runs the filtering stage of preprocessing.
            
            
            
        end
        
        function run_preprocessing(self)
            % run_preprocessing(self): Runs the preprocessing stage of the
            % analysis.
            %
            % Preprocessing turns image data (MxNxt) into a data matrix 
            % (M*Nxt).
            disp('run_preprocessing');
            
            % @todo: checks on current state
            
            set(self.h.stage_status, 'String', 'Running Pre-Processing...');
            drawnow;
            
            s = self.settings.preprocessing;
            
            if isfield(self.pp, 'im_data_mc')
                % then we have motion corrected data.
                im_data = self.pp.im_data_mc;
                im_x = self.pp.mc_x;
                im_y = self.pp.mc_y;
                im_z = 1:size(im_data, 3);
            else
                % Then use the original data.
                im_data = self.im_data;
                im_x = 1:size(im_data, 2);
                im_y = 1:size(im_data, 1);
                im_z = 1:size(im_data, 3);
            end
            
            block = ones(s.smooth_window) ./ (prod(s.smooth_window));
            
            im_conv = convn(im_data, block, 'same');
            
            self.pp.im_y = im_y(1:s.down_sample(1):length(im_y));
            self.pp.im_x = im_x(1:s.down_sample(2):length(im_x));
            self.pp.im_z = im_z(1:s.down_sample(3):length(im_z));
            
            self.pp.im_data = im_conv(self.pp.im_y, self.pp.im_x, self.pp.im_z);
            self.pp.data = reshape(self.pp.im_data, [], size(self.pp.im_data, 3));   
            
            self.pp.use_wavelets = get(self.h.wavelets_check, 'Value');
            
            if self.pp.use_wavelets
                % Then we want to do the wavelet decomposition
                disp('Computing Wavelets...');
                set(self.h.stage_status, 'String', 'Computing Wavelets...');
                drawnow;
                
                cwt_struct = cwtft(self.pp.data(1,:));
                num_scales = length(cwt_struct.scales);
                num_frames = size(self.pp.data, 2);
                
                self.pp.cwt_data = zeros(size(self.pp.data, 1), 2 * num_scales * num_frames);
                
                rcwt = real(cwt_struct.cfs);
                icwt = imag(cwt_struct.cfs);
                
                self.pp.cwt_data(1, :) = reshape([rcwt, icwt], 1, []);
                
                self.pp.cwt_struct = cwt_struct;
                
                for i = 2:size(self.pp.data, 1)
                    cwt_struct = cwtft(self.pp.data(i, :));
                    
                    rcwt = real(cwt_struct.cfs);
                    icwt = imag(cwt_struct.cfs);
                    
                    self.pp.cwt_data(i, :) = reshape([rcwt, icwt], 1, []);
                end
                
                % @todo: other things with the wavelets, like removing
                % frequencies etc. Things that can be done to cwt_data
            end

            %%% I'm going ot hack in a filter here to test
            %hd = load('icp_filter.mat');
            %data = self.pp.data - repmat(mean(self.pp.data, 2), 1, size(self.pp.data, 2));
            %self.pp.data = filter(hd.Hd, data')';            
            %%%
            
            self.stage = self.PreProcessingStage;
            
            self.update();
        end
        
        function run_pca(self)
            % run_pca(self): Runs the pca stage of the analysis.
            
            disp('run_pca');
            
            % @todo: checks on state
            if self.stage == self.InitStage
                disp('Running previous stage...');
                self.run_preprocessing();
            end
            
            set(self.h.stage_status, 'String', 'Running PCA...');
            drawnow;
            
            s = self.settings.pca;
            
            if self.pp.use_wavelets
                [scores, pcs, eigs] = princomp(self.pp.cwt_data');
            else
                [scores, pcs, eigs] = princomp(self.pp.data');
            end
            
            self.pca.scores = scores;
            self.pca.pcs = pcs;
            self.pca.eigs = eigs;
            
            
            self.pca.im_data = reshape(self.pca.scores, ...
                size(self.pp.im_data,1), size(self.pp.im_data,2), size(self.pca.scores, 2));
            
            self.pca.im_x = self.pp.im_x;
            self.pca.im_y = self.pp.im_y;
            
            self.stage = self.PcaStage;
            
            self.update();
        end
        
        function run_ica(self)
            % run_ica(self): Runs the ica stage of the analysis.
            
            disp('run_ica');
            
            % @todo: checks on state.
            if self.stage == self.InitStage || self.stage == self.PreProcessingStage
                disp('Running previous stages...');
                self.run_pca();
            end
            
            set(self.h.stage_status, 'String', 'Running ICA...');
            drawnow;
            
            s = self.settings.ica;
            
            % Make sure that the pcs requested are in range.
            which_pcs = s.which_pcs;
            which_pcs(which_pcs > size(self.pca.scores, 2)) = [];
            
            if strcmp(s.ica_func, 'CellsortICA')
                disp('Running CellsortICA');
                [ics, im_data, A, niter] = CellsortICA(self.pca.pcs(:, which_pcs)', ...
                    self.pca.im_data(:, :, which_pcs), self.pca.eigs(which_pcs), [], s.mu);
                
                self.ica.ics = ics';
                self.ica.im_data = shiftdim(im_data, 1);
                self.ica.A = A;
                self.ica.scores = reshape(self.ica.im_data, [], size(self.ica.im_data, 3));
                
            elseif strcmp(s.ica_func, 'fastica')
                disp('Running fastica');

                if s.init_guess == -1
                    % This means that the init guess should be the previous A
                    % matrix.
                    if isfield(self.ica, 'A') && size(self.ica.A, 1) == length(which_pcs)
                        % Then the last A matrix is valid.
                        [icasig, A, W] = fastica(self.pca.scores(:, which_pcs)', 'initGuess', self.ica.A);
                    else
                        % Then the last A matrix is invalid.
                        disp('Could not use previous A as initial guess.');                        
                        [icasig, A, W] = fastica(self.pca.scores(:, which_pcs)');
                    end
                elseif s.init_guess == 0
                    % Then we just use the normal random guess.                        
                    %%%
                    %[icasig, A, W] = fastica(self.pca.scores(:, which_pcs)');
                    %%% Messing around with the other fastica params
                    disp('fastica params');
                    [icasig, A, W] = fastica(self.pca.scores(:, which_pcs)', 'g', 'tanh', 'a1', 3, 'stabilization', 'on');
                    %%%
                else
                    % Then s.init_guess should be the A matrix itself.
                    if size(s.init_guess, 1) == length(which_pcs)
                        disp('Using Init Guess');
                        [icasig, A, W] = fastica(self.pca.scores(:, which_pcs)', 'initGuess', s.init_guess);
                    else
                        disp('Initial Guess incorrectly formatted.');
                        [icasig, A, W] = fastica(self.pca.scores(:, which_pcs)');
                    end
                end
                
                self.ica.scores = icasig';

                if s.positive_skew
                    % Then we want the component scores to all skew positively.
                    is_neg_skew = skewness(self.ica.scores) < 0;
                    self.ica.scores(:, is_neg_skew) = -self.ica.scores(:, is_neg_skew);
                    A(:, is_neg_skew) = -A(:, is_neg_skew);
                end
                
                self.ica.A = A;
                self.ica.W = W;
                ics = self.pca.pcs(:, which_pcs) * A;
                
                if self.pp.use_wavelets
                    % Reconstruct the original data from wavelets
                    cwt_struct = self.pp.cwt_struct;
                    num_frames = size(self.pp.data, 2);
                    wavelet_rec = zeros(num_frames, size(ics, 2));
                    for i = 1:size(ics, 2)
                        ic_cfs = reshape(ics(:,i), length(cwt_struct.scales), []);
                        cwt_struct.cfs = complex(ic_cfs(:, 1:num_frames), ic_cfs(:, (num_frames+1):end)); 
                        
                        wavelet_rec(:, i) = icwtft(cwt_struct);
                    end
                    self.ica.ics = wavelet_rec;
                else
                    self.ica.ics = ics;
                end
               
                
                
                self.ica.im_data = reshape(self.ica.scores, ...
                    size(self.pp.im_data,1), size(self.pp.im_data,2), size(self.ica.scores, 2));
            elseif strcmp(s.ica_func, 'TreeICA')
                
            elseif strcmp(s.ica_func, 'RadICAl')
                
            elseif strcmp(s.ica_func, 'KernelICA')
                disp('Running Kernel ICA');
                tic;
                W = kernel_ica(self.pca.scores(:, which_pcs)');
                
                self.ica.scores = (W * self.pca.scores(:, which_pcs)')';
                self.ica.ics = self.pca.pcs(:, which_pcs) / W;
                
                self.ica.im_data = reshape(self.ica.scores, ...
                    size(self.pp.im_data,1), size(self.pp.im_data, 2), size(self.ica.scores, 2));
                
                toc;
            elseif strcmp(s.ica_func, 'imageica')
                
            elseif strcmp(s.ica_func, 'icasso')
                
            elseif strcmp(s.ica_func, 'fourierica')
                [S_Ft, A, W] = fourierica(self.pca.scores(:, which_pcs)', length(which_pcs), 50, 0, 20);
                
                self.ica.scores = (W * self.pca.scores(:, which_pcs)')';
                self.ica.ics = self.pca.pcs(:, which_pcs) * A;
                
                self.ica.im_data = reshape(self.ica.scores, ...
                    size(self.pp.im_data,1), size(self.pp.im_data, 2), size(self.ica.scores, 2));
            else
                error('Incorrect ICA function');
            end
            self.ica.im_x = self.pca.im_x;
            self.ica.im_y = self.pca.im_y;
            
            self.stage = self.IcaStage;
            
            self.update();
        end
        
        function run_segmentation(self)
            % run_segmentation(self): Runs the segmentation algorithm on
            % the ICs.
            %
            % This looks for ic components that are spatially localized.
            
            if self.stage < self.IcaStage
                disp('Run ICA first');
                return;
            end
            
            set(self.h.stage_status, 'String', 'Running Segmentation...');
            drawnow;
            
            [segment_masks, segment_info] = segment_ics(self.ica.im_data, self.ica.im_x, self.ica.im_y);
            
            self.post.im_data = segment_masks;
            self.post.im_x = self.ica.im_x;
            self.post.im_y = self.ica.im_y;
            self.post.segment_info = segment_info;
            
            self.stage = self.PostProcessingStage;
            self.update();
        end
        
        function rois = calc_rois(self)
            % calc_rois(self): calculates the ROIs from the independent
            % components.
            
            if self.stage < self.IcaStage
                disp('Run ICA first');
                return;
            end
            
            
            rois = self.settings.post.calc_rois_func(self.ica.im_data, self.ica.im_x, self.ica.im_y);
            
            self.create_new_roi_set('IC_Generated_ROIs', true);
            self.rois{self.gui.current_roi_set} = rois;
            
            self.h.roi_editor.set_rois(self.rois{self.gui.current_roi_set}, false);
            
            self.update();
            
        end
        
        function estimate_clusters(self)
            disp('estimate_clusters');
            
            if self.stage < self.IcaStage
                disp('Run ICA first');
                return;
            end
            
            % We need rois specifically for clustering.
            if self.gui.is_component_set(self.gui.current_roi_set)
                self.cluster.rois = self.rois{self.gui.current_roi_set};
            else
                self.cluster.rois = self.calc_rois();
            end
            
            s = self.settings.cluster;
            
            [feature_matrix, feature_weights] = s.feature_matrix_func(self.ica.ics, self.ica.im_data);
            
            self.cluster.feature_matrix = feature_matrix;
            self.cluster.feature_weights = feature_weights;
            
            self.cluster.similarity_matrix = s.calc_sim_matrix_func(self.cluster.feature_weights, self.cluster.feature_matrix);
            
            [self.cluster.im_data, top3] = s.visualization_func(self.cluster.similarity_matrix, self.ica.im_data);
            
            self.cluster.im_x = self.ica.im_x;
            self.cluster.im_y = self.ica.im_y;
            
            self.cluster.top3 = top3;
            
            %%% Not sure how/when to set this
            self.cluster.cluster_matrix = zeros(size(self.cluster.similarity_matrix));
            % THis is just a hack for testing
            self.cluster.cluster_matrix(1, 5) = 1;
            self.cluster.cluster_matrix(2, [10, 20]) = 1;
            %%%
            
            self.stage = self.PostProcessingStage;
            self.update();
        end
        
        function load_settings(self, filename)
            % load_settings(self, filename): loads settings from a saved
            % file.
            
            disp('load_settings');
        end
        
        function save_settings(self, filename)
            % save_settings(self, filename): saves the settings to a file.
            disp('save_settings');
        end
        
        function edit_preprocessing_settings(self)
            % edit_preprocessing_settings(self): shows the edit
            % preprocessing settings dialog.
            
            self.h.pp_dialog = dialog('Name', 'Edit Preprocessing Settings', 'Units', 'pixels');
            
            self.h.pp_main_vbox = uiextras.VBox('Parent', self.h.pp_dialog, ...
                'Spacing', self.gui.MARGIN, 'Padding', self.gui.MARGIN);
            self.h.pp_settings_grid = uiextras.Grid('Parent', self.h.pp_main_vbox, 'Spacing', self.gui.MARGIN);
            self.h.pp_button_hbox = uiextras.HButtonBox('Parent', self.h.pp_main_vbox, 'Spacing', self.gui.MARGIN);
            
            % The order of these matters, as they are added to the grid
            % based on the order.
            % First column
            uiextras.Empty('Parent', self.h.pp_settings_grid);
            self.h.pp_smooth_label = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'text', 'String', 'Smooth Window:', 'TooltipString', ...
                'Sets dimensions of smooth window - MxNxT', ...
                'HorizontalAlignment', 'left');
            self.h.pp_down_label = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'text', 'String', 'Down Sampling:', 'TooltipString', ...
                'Sets the down sampling in each dimensoin - MxNxT', ...
                'HorizontalAlignment', 'left');
            % Second column
            self.h.pp_M_label = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'text', 'String', 'M', 'TooltipString', ...
                'M corresponds to the rows of the image (vertical)');
            self.h.pp_smooth_M = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'edit', 'String', self.settings.preprocessing.smooth_window(1));
            self.h.pp_down_M = uicontrol('Parent', self.h.pp_settings_grid,...
                'Style', 'edit', 'String', self.settings.preprocessing.down_sample(1));
            
            % Third column
            self.h.pp_N_label = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'text', 'String', 'N', 'TooltipString', ...
                'N corresponds to the columns of th eimage (horizontal)');
            self.h.pp_smooth_N = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'edit', 'String', self.settings.preprocessing.smooth_window(2));
            self.h.pp_down_N = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'edit', 'String', self.settings.preprocessing.down_sample(2));
            
            % Fourth column
            self.h.pp_T_label = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'text', 'String', 'T', 'TooltipString', ...
                'T corresponds to the frames (depth)');
            self.h.pp_smooth_T = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'edit', 'String', self.settings.preprocessing.smooth_window(3));
            self.h.pp_down_T = uicontrol('Parent', self.h.pp_settings_grid, ...
                'Style', 'edit', 'String', self.settings.preprocessing.down_sample(3));
            
            
            self.h.pp_ok_button = uicontrol('Parent', self.h.pp_button_hbox, ...
                'Style', 'pushbutton', 'String', 'OK', 'Callback', @self.pp_settings_ok_cb);
            self.h.pp_cancel_button = uicontrol('Parent', self.h.pp_button_hbox, ...
                'Style', 'pushbutton', 'String', 'Cancel', 'Callback', @self.pp_settings_cancel_cb);
            
            set(self.h.pp_settings_grid, 'ColumnSizes', [-1, 30, 30, 30], 'RowSizes', [30, 30, 30]);
            set(self.h.pp_main_vbox, 'Sizes', [-1, 30]);
            
            fh_pos = get(self.h.fh, 'Position');
            set(self.h.pp_dialog, 'Position', [fh_pos(1) + 100, fh_pos(2) + 100, 300, 200]);
        end
        
        function set_pp_smooth_window(self, smM, smN, smT)
            % set_smooth_window: sets the size of the smoothing convolution
            % window.
            
            if nargin == 2
                % Then the 3 values should all be in smM
                assert(length(smM) == 3);
                smN = smM(2);
                smT = smM(3);
                smM = smM(1);
            elseif nargin < 2 || isempty(smM)
                % Default will be to leave unchanged.
                smM = NaN;
            elseif nargin < 3 || isempty(smN)
                smN = NaN;
            elseif nargin < 4 || isempty(smT)
                smT = NaN;
            end
            
            % Get the old values
            sm_old = self.settings.preprocessing.smooth_window;
            
            % Set the new values
            self.settings.preprocessing.smooth_window = [smM, smN, smT];
            
            % Check for bad values
            nan_vals = isnan(self.settings.preprocessing.smooth_window) ...
                | self.settings.preprocessing.smooth_window < 0;
            if any(nan_vals)
                % Then some of the values are bad/missing, use the old
                % values.
                disp('Some smooth window values bad/missing.');
                self.settings.preprocessing.smooth_window(nan_vals) = sm_old(nan_vals);
            end
        end
        
        function set_pp_down_sample(self, dsM, dsN, dsT)
            % set_down_sample: sets the size of the down sampling
            
            if nargin == 2
                % Then the 3 values should all be in smM
                assert(length(dsM) == 3);
                dsN = dsM(2);
                dsT = dsM(3);
                dsM = dsM(1);
            elseif nargin < 2 || isempty(dsM)
                % Default will be to leave unchanged.
                dsM = NaN;
            elseif nargin < 3 || isempty(dsN)
                dsN = NaN;
            elseif nargin < 4 || isempty(dsT)
                dsT = NaN;
            end
            
            % Get the old values
            ds_old = self.settings.preprocessing.down_sample;
            
            % Set the new values
            self.settings.preprocessing.down_sample = [dsM, dsN, dsT];
            
            % Check for bad values
            nan_vals = isnan(self.settings.preprocessing.down_sample) ...
                | self.settings.preprocessing.down_sample < 0;
            if any(nan_vals)
                % Then some of the values are bad/missing, use the old
                % values.
                disp('Some down sample values bad/missing.');
                self.settings.preprocessing.down_sample(nan_vals) = ds_old(nan_vals);
            end
        end
        
        function edit_pca_settings(self)
            
        end
        
        function edit_ica_settings(self)
            % edit_ica_settings(self): shows the edit ica settings dialog.
            
            self.h.ica_dialog = dialog('Name', 'Edit Preprocessing Settings', 'Units', 'pixels');
            
        end
        
        function set_data(self, im)
            % set_data(self, im): sets the image data.
            %
            % @param: im image data (MxNxt)
            disp('set_data');
            
            % @todo: checks on inputs
            
            self.im_data = im;
            self.pp = [];
            self.pca = [];
            self.ica = [];
            
            self.stage = self.InitStage;
            
            self.update();
        end
        
        function set.Parent(self, val)
            disp('Setting RoiEditor Parent');
            set(self.h.panel, 'Parent', double(val))
        end
        
        function component = get_current_component(self)
            % Returns the currently selected component
            component = [];
            f = self.h.roi_editor.current_frame;

            switch self.stage
                case self.InitStage
                    return
                case self.PreProcessingStage
                    return
                case self.PcaStage
                    component = self.pca.pcs(:, f);
                    return
                case self.IcaStage
                    component = self.ica.ics(:, f);
            end
        end
        
        function lock_current_component(self)
            c = self.get_current_component();
            
            if ~isempty(c)
                self.gui.locked_components{end+1} = c;
            else
                disp('No component to lock.');
            end
            
            %self.update();
        end
        
        function clear_locked_components(self)
            self.gui.locked_components = {};
            self.update();
        end
        
        function show_best_selected_components(self, rois)
            % This plots the top components with an roi 
            
            disp('show_best_selected_components');
            
            if isempty(rois)
                % Then nothing to plot.
                return
            end
            
            im_data = [];
            comp = [];
            im_x = [];
            im_y = [];
            
            if self.gui.display == 4 && self.stage == self.IcaStage
                % Then view the ICs
                im_data = self.ica.im_data;
                im_x = self.ica.im_x;
                im_y = self.ica.im_y;
                comp = self.ica.ics;
            elseif self.gui.display == 3 && self.stage >= self.PcaStage
                % Then view the PCs
                im_data = self.pca.im_data;
                im_x = self.pca.im_x;
                im_y = self.pca.im_y;
                comp = self.pca.pcs;
            elseif self.stage == self.PcaStage
                % Then view the PCs
                im_data = self.pca.im_data;
                im_x = self.pca.im_x;
                im_y = self.pca.im_y;
                comp = self.pca.pcs;
            elseif self.stage == self.IcaStage
                % Then view the ICs
                im_data = self.ica.im_data;
                im_x = self.ica.im_x;
                im_y = self.ica.im_y;
                comp = self.ica.ics;
            else
                return;
            end

            s = self.settings.visualize;
            
            top_idxs = [];
            top_vals = [];
            residuals = [];
            for i = 1:size(rois, 1)
                % Lets just start with the center of the roi.
                cx = interp1(im_x, 1:size(im_data, 2), rois(i, 1), 'nearest');
                cy = interp1(im_y, 1:size(im_data, 1), rois(i, 2), 'nearest');
                
                % @todo: check if in the right range
                
                % Now get the 3 largest scores at the point
                point_data = squeeze(im_data(cy, cx, :));
                
%                 if s.best_selection_type == 0
%                     % Then use the data exactly
%                     
%                 elseif s.best_selection_type == 1
%                     % Use the absolute value
%                     point_data = abs(point_data);
%                 elseif s.best_selection_type == 2
%                     % 
%                 end
                
                [vals, sidx] = sort(point_data, 1, 'descend');
                
                n_best = min(s.num_best_components, length(sidx));
                top_idxs = sidx(1:n_best);
                top_vals = vals(1:n_best);
                
                if length(sidx > n_best)
                    residual_idxs = sidx((n_best+1):end);
                    residuals(i,:) = sum(comp(:, residual_idxs) .* repmat(vals((n_best+1):end)', size(comp, 1), 1), 2);
                end
            end
            
            [top_idxs, ia, ic] = unique(top_idxs);
            top_vals = top_vals(ia);
            
            cla(self.h.component_axes);
            legend(self.h.component_axes, 'off');
            hold(self.h.component_axes, 'all');
            
%             if isfield(self.h, 'component_labels')
%                 delete(self.h.component_labels);
%             end
            
            self.h.component_labels = [];
            
            for i = 1:length(top_idxs)
                c = self.gui.locked_colors(mod(i-1, length(self.gui.locked_colors))+1, :); 
                
                component = comp(:, top_idxs(i));
                
                if ~s.norm_best_components
                    % Then multiply the component by its value.
                    component = component * top_vals(i);
                end
                
                plot(self.h.component_axes, component, 'Color', c);
                text_x = 0.95 - (length(top_idxs) + size(residuals, 1) - i) * 0.05;
                self.h.component_labels(i) = text(text_x, 0.9, num2str(top_idxs(i)), ...
                    'Parent', self.h.component_axes, 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'UserData', top_idxs(i), ...
                    'Color', c, 'ButtonDownFcn', @self.component_label_cb);
            end

            if s.show_residual
                for i = 1:size(residuals, 1)
                    c = self.gui.locked_colors(mod(i+length(top_idxs)-1, length(self.gui.locked_colors))+1, :);

                    plot(self.h.component_axes, residuals(i,:), ':', 'Color', c);
                    text_x = 0.95 - (size(residuals, 1) - i) * 0.05;
                    self.h.component_labels(i) = text(text_x, 0.9, 'res', ...
                        'Parent', self.h.component_axes, 'Units', 'normalized', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', c, 'ButtonDownFcn', @self.component_label_cb);
                end
            end
            
            %self.h.lh = legend(self.h.component_axes, 'show');
        end
        
        function create_new_roi_set(self, name, is_component_set)
            % create_new_roi_set(self, name): creates a new roi set.
            
            disp('create_new_roi_set');
            
            if nargin < 3 || isempty(is_component_set)
                is_component_set = false;
            end
            
            if nargin < 2 || isempty(name)
                % Then prompt for a name for the roi set.
                name = inputdlg('Enter name for ROI set:');
                name = name{1};
            end
            
            if isempty(name)
                % Then at this point it was canceled
                return;
            end
            
            % @todo: other checks on name
            
            % Add a new empty set of rois.
            self.rois{end+1} = [];
            self.gui.current_roi_set = length(self.rois);
            self.gui.roi_names{length(self.rois)} = name;
            self.gui.is_component_set(length(self.rois)) = is_component_set;
            
            self.h.roi_editor.set_rois(self.rois{self.gui.current_roi_set});
            
            self.update();
        end
        
        function delete_roi_set(self, idx)
            % delete_roi_set(self, idx): deletes an roi set.
            
            disp('delete_roi_set');
            
            if nargin < 2
                % Then just delete the currently selected set(?)
                idx = self.gui.current_roi_set;
            end
            
            % prompt?
            
            % Make sure idx is in range
            if idx < 1 || idx > length(self.rois)
                % umm... do nothing?
                disp('Requested ROI set not found. Nothing deleted.');
                return;
            end
            
            self.rois(idx) = [];
            self.gui.roi_names(idx) = [];
            self.gui.is_component_set(idx) = [];
            
            if isempty(self.rois)
                self.build_roi_listbox();
            end
            
            if self.gui.current_roi_set > length(self.rois)
                self.gui.current_roi_set = length(self.rois);
            end
            
            self.h.roi_editor.set_rois(self.rois{self.gui.current_roi_set}, false);
            self.update();
        end
        
        function select_roi_set(self, idx)
            % select_roi_set(self, idx): sets the currently selected roi
            % set.
            
            disp('select_roi_set');
            
            if nargin < 2 || isempty(idx)
                idx = 1;
            end
            
            % Check for empty
            if isempty(self.rois)
                self.build_roi_listbox();
            end
            
            % Make sure idx is in range.
            if idx < 1 || idx > length(self.rois)
                % Then idx is out of range.
                idx = length(self.rois);
            end
            
            if idx == self.gui.current_roi_set
                % Then nothing needs to be done.
                return;
            end
            
            self.gui.current_roi_set = idx;
            
            self.h.roi_editor.set_rois(self.rois{self.gui.current_roi_set}, false);
            
            self.update();
        end
        
        function set_ica_spatial_guess(self, masks, which_pcs)
            % set_ica_spatial_guess(self, masks): Sets an initial guess for
            % the ICs given a series of spatial masks.
            %
            % This adjusts the settings.ica.which_pcs to match the number
            % of masks.
            %
            % @param: masks mxnxG array of masks indicating the location of
            % the ics. G is the number of masks.
            % @param: which_pcs the pcs to use for the guess, this must
            % have a length of G. Default is to use 1:G pcs
            
            if nargin < 3 || isempty(which_pcs)
                which_pcs = 1:size(masks, 3);
            end
            
            % Must have already run pca at least.
            if self.stage == self.InitStage || self.stage == self.PreProcessingStage
                disp('Running previous stages...');
                self.run_pca();
            end
            
            pc_scores = self.pca.scores(:, which_pcs)';
            
            % reformat the masks into scores form
            ic_scores = reshape(masks, [], size(masks, 3))';
            
            A_guess = pc_scores / ic_scores;
            
            self.settings.ica.init_guess = A_guess;
            self.settings.ica.which_pcs = which_pcs;
        end
        
        function save_rois(self, filename)
            % save_rois(self, filename): saves the rois to a mat file.
            % 
            % If filename is not given then user will be prompted
            
            if nargin < 2 || isempty(filename)
                % Prompt for filename
                [file, path] = uiputfile({'*.mat'}, 'Save ROIs');
                if file == 0
                    % Dialog was canceled
                    return;
                end
                filename = [path file];
            end
            rois = self.rois;
            names = self.gui.roi_names;
            
            save(filename, 'rois', 'names');
        end
        
        function load_rois(self, filename)
            % load_rois(self, filename): loads rois from a mat file.
            
            if nargin < 2 || isempty(filename)
                [file, path] = uigetfile({'*.mat'}, 'Load ROIs');
                
                if file == 0
                    % Dialog was canceled
                    return;
                end
                filename = [path file];
            end
            
            roi_s = load(filename);
            
            for i = 1:length(roi_s.rois)
                self.create_new_roi_set(roi_s.names{i});
                self.rois{self.gui.current_roi_set} = roi_s.rois{i};
            end
            
            self.h.roi_editor.set_rois(self.rois{self.gui.current_roi_set}, false);
            
            self.update();
        end
        
        function set_cluster(self, ic_idxs)
            % set_cluster(self, ic_idxs): Sets the ics given by ic_idxs as
            % a cluster.
            
            disp('set_cluster');
            
            % @todo: checks on state, cluster
            if nargin < 2 || isempty(ic_idxs)
                % Then the ic_idxs should be the rois that are selected.
                ic_idxs = self.h.roi_editor.selected_rois;
            end
            % @todo: more checks
            
            % Basically, just set the cluster matrix to 1 at the ic_idxs.
            self.cluster.cluster_matrix(ic_idxs, ic_idxs) = 1;
            
            self.update();
        end
        
        function uncluster(self, idxs)
            % uncluster(self, idxs): unclusters the given idxs.
            disp('uncluster');
            
            if nargin < 2 || isempty(idxs)
                % Then idxs are the currently selected rois
                idxs = self.h.roi_editor.selected_rois;
            end
            
            self.cluster.cluster_matrix(idxs, idxs) = 0;
            
            self.update();
        end
        
        function clear_clusters(self)
            % clear_clusters(self): removes all clusters
        end
                
        %%%%%%%%%% Utility Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function default_settings(self)
            % default_settings(self): sets the settings to the default.
            
            %self.settings.preprocessing.block_size = 6;
            self.settings.preprocessing.smooth_window = [6, 6, 1];
            %self.settings.preprocessing.step_size = 2;
            self.settings.preprocessing.down_sample = [2, 2, 1];
            
            self.settings.preprocessing.motion_correct_func = @align_im_stack;
            
            self.settings.preprocessing.default_filter_file = 'icp_default_filters.mat';
            self.settings.preprocessing.filter_file = self.settings.preprocessing.default_filter_file;
            
            self.settings.pca.mean_center = 0;
            
            self.settings.ica.which_pcs = 1:150;
            self.settings.ica.positive_skew = true;
            self.settings.ica.init_guess = -1;
            self.settings.ica.ica_func = 'fastica';
            self.settings.ica.mu = 0.2; % This is a parameter for Cellsort ICA
            
            self.settings.post.calc_rois_func = @calc_rois_from_components;
            self.settings.post.segment_ics_func = @segment_ics;
            
            self.settings.cluster.feature_matrix_func = @make_feature_matrix;
            self.settings.cluster.feature_weights = 1;
            self.settings.cluster.visualization_func = @viz_best3;
            self.settings.cluster.calc_sim_matrix_func = @(b, x) sum(repmat(b, [size(x,1), size(x,2), 1]) .* x, 3);
            
            
            self.settings.visualize.norm_best_components = false;
            self.settings.visualize.num_best_components = 3;
            self.settings.visualize.show_residual = true;
            self.settings.visualize.best_selection_type = 1;            
            
            
        end
        
        function plot_locked_components(self)
            hold(self.h.component_axes, 'on');
            for i = 1:length(self.gui.locked_components)
                c = self.gui.locked_colors(mod(i-1, length(self.gui.locked_colors))+1, :); 
                plot(self.h.component_axes, self.gui.locked_components{i}, 'Color', c);
            end
        end
        
        function plot_3d(self)
            % plot_3d(self): plots the currently selected components in
            % 3-space.
            disp('plot_3d');

            if self.gui.display == 4 && self.stage == self.IcaStage
                % Then view the ICs
                sc = self.ica.scores;
                self.build_3d_gui_components(self.IcaStage)
            elseif self.gui.display == 3 && self.stage >= self.PcaStage
                % Then view the PCs
                sc = self.pca.scores;
                self.build_3d_gui_components(self.PcaStage);

            elseif self.stage == self.PcaStage
                % Then view the PCs
                sc = self.pca.scores;
                self.build_3d_gui_components(self.PcaStage);
            elseif self.stage == self.IcaStage
                % Then view the ICs
                sc = self.ica.scores;
                self.build_3d_gui_components(self.IcaStage);
            else
                cla(self.h.pc3_axes);
                self.build_3d_gui_components(self.stage);
                return;
            end
            
            if isempty(sc)
                % Then there's nothing to plot
                return;
            end
            
            % Make sure everything is in range
            if self.gui.x_pc < 1 || self.gui.x_pc > size(sc, 2)
                % Just set to 1 for now.
                self.gui.x_pc = 1;
            end
            if self.gui.y_pc < 1 || self.gui.y_pc > size(sc, 2)
                self.gui.y_pc = 1;
            end
            if self.gui.z_pc < 1 || self.gui.z_pc > size(sc, 2)
                self.gui.z_pc = 1;
            end
            
            xv = sc(:, self.gui.x_pc);
            yv = sc(:, self.gui.y_pc);
            zv = sc(:, self.gui.z_pc);
            
            [az, el] = view(self.h.pc3_axes);
            plot3(self.h.pc3_axes, xv, yv, zv, '.k');
            axis(self.h.pc3_axes, 'vis3d');
            view(self.h.pc3_axes, az, el);
            
            popup_str = get(self.h.x_popup, 'String');
            xlabel(self.h.pc3_axes, popup_str{self.gui.x_pc});
            ylabel(self.h.pc3_axes, popup_str{self.gui.y_pc});
            zlabel(self.h.pc3_axes, popup_str{self.gui.z_pc});
        end
        
        function build_3d_gui_components(self, stage)
            % build_3d_gui_components(self): builds the 3d axis gui
            
            switch stage
                case {self.InitStage, self.PreProcessingStage}
                    % There's nothing to do in this case
                    set(self.h.x_popup, 'String', {'None'}, 'Value', 1);
                    set(self.h.y_popup, 'String', {'None'}, 'Value', 1);
                    set(self.h.z_popup, 'String', {'None'}, 'Value', 1);
                    return;
                case self.PcaStage
                    popup_str = cellfun(@(n) ['PC ' num2str(n)], num2cell(1:size(self.pca.scores, 2)), ...
                        'UniformOutput', 0);
                    set(self.h.x_popup, 'String', popup_str);
                    set(self.h.y_popup, 'String', popup_str);
                    set(self.h.z_popup, 'String', popup_str);                    
                case self.IcaStage
                    popup_str = cellfun(@(n) ['IC ' num2str(n)], num2cell(1:size(self.ica.scores, 2)), ...
                        'UniformOutput', 0);
                    set(self.h.x_popup, 'String', popup_str);
                    set(self.h.y_popup, 'String', popup_str);
                    set(self.h.z_popup, 'String', popup_str);                    
            end
            
            set(self.h.x_popup, 'Value', self.gui.x_pc);
            set(self.h.y_popup, 'Value', self.gui.y_pc);
            set(self.h.z_popup, 'Value', self.gui.z_pc);
        end
        
        function build_roi_listbox(self)
            % build_roi_listbox(self): creates the roi_listbox.
            
            if isempty(self.rois)
                % Then there is nothing in the rois, so create an empty
                % list.
                self.rois{1} = [];
                self.gui.current_roi_set = 1;
                self.gui.roi_names = {};
                self.gui.roi_names{1} = 'ROI_set_1';
                self.gui.is_component_set = false;
            end
            
            % Check ranges
            assert(length(self.rois) == length(self.gui.roi_names));
            if self.gui.current_roi_set < 1 || self.gui.current_roi_set > length(self.rois)
                self.gui.current_roi_set = length(self.rois);
            end
            
            set(self.h.roi_listbox, 'String', self.gui.roi_names, 'Value', self.gui.current_roi_set);
        end
        
        function update_current_roi_set(self)
            % update_current_roi_set(self): synchronizes the current roi
            % set with the roi_editor.
            
            disp('update_current_roi_set');
            
            self.rois{self.gui.current_roi_set} = self.h.roi_editor.rois.xyrra(:,:,1);
            self.update();
        end
        
        function draw_clustered_rois(self)
            % Ok
            disp('draw_clustered_rois');
            
            roi_idx = self.h.roi_editor.selected_rois;
            roi_h = self.h.roi_editor.h.rois;
            
            set(roi_h, 'Color', [0.5, 0.5, 0.5], 'Marker', 'none');
            
            set(roi_h(roi_idx), 'Color', 'w');
            
            % First see if there are other clustered rois
            cluster_idx = [];
            if ~isempty(roi_idx)
                cluster_idx = find(self.cluster.cluster_matrix(roi_idx(1), :));
            end
            
            if ~isempty(cluster_idx)
                % Then highlight the clustered rois
                set(roi_h(cluster_idx), 'Color', [0.9, 0.9, 0.9]);                
            end
            
            % Draw the estimates
            
            f = self.h.roi_editor.current_frame;
            
            set(roi_h(f), 'Marker', 's', 'MarkerSize', 1, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'w');
            set(roi_h(self.cluster.top3(f, 1)), 'Marker', 's', 'MarkerSize', 1, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            set(roi_h(self.cluster.top3(f, 2)), 'Marker', 's', 'MarkerSize', 1, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
            set(roi_h(self.cluster.top3(f, 3)), 'Marker', 's', 'MarkerSize', 1, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');

            
        end
    end
    
end