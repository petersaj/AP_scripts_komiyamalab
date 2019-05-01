classdef AxesEventNotifier < handle
    % class AxesEventNotifier: sends out notifications during axes events.
    
    properties
        ah; % The axes handles associated with this notifier.
        source_h; % The handle of the axes where the event occured
        button_down; % Which mouse button is down.
        last_eventdata; % Keeps the data given by eventdata from the last event.
    end
    
    events
        ButtonDown;
    end
    
    methods 
        function self = AxesEventNotifier(ah)
            % Constructor
            self.ah = [];
            % First check if there is already a AxesEventNotifier for
            % this axes.
            user_data = get(ah(1), 'UserData');
            if ~isempty(user_data)
                if isa(user_data, 'AxesEventNotifier')
                    % Then the axes already has an AxesEventNotifier.
                    self = user_data;
                    return;
                else
                    % Then there is something else in user_data. Issue a
                    % warning.
                    warning('AxesEventNotifier:deletingUserData', ...
                        'Axes UserData is not empty and will be replaced.');
                end
            end
            
            set(ah, 'UserData', self);
            self.add_handle(ah);
        end
        
        %%%%%%%%%% Callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function buttondown_cb(self, source_h, eventdata)
            % Called when a button is pressed in an axes.            
            self.source_h = source_h;
            self.last_eventdata = eventdata;
            self.notify('ButtonDown');
        end
        
        %%%%%%%%%% Main Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function add_handle(self, ah)
            set(ah, 'ButtonDownFcn', @self.buttondown_cb);
            self.ah = [self.ah ah];
        end
    end
end