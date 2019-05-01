function AP_labjack_lever_single_threshold
% Run program to detect lever presses via labjack/lever
% Single threshold, defined by distance and time in GUI
% Requires touch sensor for lever - presses only triggered if hand on
% handle


clear all

%% Initialize LabJack

% For some reason it usually takes 2 tries to initalize labjack, this just
% does that automatically
labjack_ready = false;
while ~labjack_ready
    try
        ljasm = NET.addAssembly('LJUDDotNet'); %Make the UD .NET assembly visible in MATLAB
        ljudObj = LabJack.LabJackUD.LJUD;
        
        %Used for casting a value to a CHANNEL enum
        chanType = LabJack.LabJackUD.CHANNEL.LOCALID.GetType;
        
        %Open the first found LabJack U3.
        [ljerror, ljhandle] = ljudObj.OpenLabJack(LabJack.LabJackUD.DEVICE.U3, LabJack.LabJackUD.CONNECTION.USB, '0', true, 0);
        
        %Start by using the pin_configuration_reset IOType so that all
        %pin assignments are in the factory default condition.
        chanObj = System.Enum.ToObject(chanType, 0); %channel = 0
        ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PIN_CONFIGURATION_RESET, chanObj, 0, 0);
        % if that throws an error, possible alt code w/o setting chanObj workaround
        %ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PIN_CONFIGURATION_RESET, 0, 0, 0);
        
        % Open FIO4-5 as AI ports (lever)
        ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_ANALOG_ENABLE_BIT,4,1,0,0);
        ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_ANALOG_ENABLE_BIT,5,1,0,0);
        % FIO7 should be set to digital by default (touch state)
        
        labjack_ready = true;
    catch me
    end
end

% Execute all queued commands
ljudObj.GoOne(ljhandle);

%Request a reading from FIO5-FIO4 (lever)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_AIN_DIFF, 5, 0, 4, 0);

% Request a reading from FIO7 (touch state)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 7, 0, 0, 0);

% Set the initial DAC0 output to 0 (trigger output)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DAC, 0, 0, 0, 0);

ioDummy = LabJack.LabJackUD.IO;
chanDummy = LabJack.LabJackUD.CHANNEL;


%% Set up buffer, start timers

% Set number of samples to plot
plot_samples = 1000;

% Kill any lingering timers (shouldn't be anything)
curr_timers = timerfindall;
if ~isempty(curr_timers)
    stop(curr_timers);
    delete(curr_timers);
end

% Define acquisition rate (in Hz)
% This can be very sensitive: if it's trying to go too fast it can get hung
% up in the timer function and not display. Can depend on computer,
% currently running processes, etc.
getLever_period = 1/50;
plotLever_period = 1/50;

% Define the acquisition/plotting timers
getLever_timer = timer('timerFcn',@getLever,'period',getLever_period,'executionMode','fixedRate');
plotLever_timer = timer('timerFcn',@plotLever,'period',plotLever_period,'executionMode','fixedRate');

% Create figure for plotting
lever_figure = figure('Position',[165    68   956   336]);

% Create plotting axes
strip_chart = subplot(4,2,[1 3]); hold on;
xlim([0 plot_samples]);
title('Lever voltage')

output_plot = subplot(4,2,5); hold on;
xlim([0 plot_samples]);
ylim([-1 6]);
title('Output voltage')

rate_plot = subplot(4,2,7); hold on;
xlim([0 plot_samples]);
title('Sampling rate offset')

% Set default values
rest_threshold_default = -2.2; % V
press_threshold_default = -2.4; % V
threshold_time_default = 200; % ms

% Get, group current trajectory information
lever_threshold.rest_threshold = rest_threshold_default;
lever_threshold.press_threshold = press_threshold_default;
lever_threshold.threshold_time = threshold_time_default;

% Set buffer sizes for analysis, initialize
% lever to keep in buffer is one more than possible to detect given time
analyze_samples = ceil(threshold_time_default)/(getLever_period*1000) + 1;
lever_analyze_buffer = nan(analyze_samples,1);
lever_time_analyze_buffer = nan(analyze_samples,1);
lever_plot_buffer = nan(plot_samples,1);
lever_rate = nan(plot_samples,1);
output_buffer = nan(plot_samples,1);

% Create an output buffer (because one sample doesn't trigger)
output_samples = 5;
output_signal = zeros(output_samples,1);

% Set initial touch state
touch_state = NaN;

% Set last success time and refractory period (in seconds)
last_success_time = clock;
refractory_period = 1;

% Create plot handles for lever and sampling rate
lever_plot_h = plot(strip_chart,lever_plot_buffer,'k');
output_plot_h = plot(output_plot,output_buffer,'k');
rate_plot_h = plot(rate_plot,lever_rate,'k');

% Create start/stop buttons
start_h = uicontrol('Style','pushbutton', ...
    'Position',[545 15 100 70],'String','Start','FontSize',30, ...
    'Callback',@start_callback);

stop_h = uicontrol('Style','pushbutton', ...
    'Position',[660 15 100 70],'String','Stop','FontSize',30, ...
    'Callback',@stop_callback);

% Create editable fields for lever threshold/floor/auto-narrow
figure_color = get(lever_figure,'Color');

uicontrol(lever_figure,'Style','text','String','Rest threshold', ...
    'BackgroundColor',figure_color,'Position',[545 120 70 30]);
edit_lever_restthresh_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(rest_threshold_default),...
    'Position',[545 100 70 30],'Callback',@updateThreshold_manual);

uicontrol(lever_figure,'Style','text','String','Press threshold', ...
    'BackgroundColor',figure_color,'Position',[675 130 70 30]);
edit_lever_pressthresh_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(press_threshold_default),...
    'Position',[675 100 70 30],'Callback',@updateThreshold_manual);

uicontrol(lever_figure,'Style','text','String','Time to threshold', ...
    'BackgroundColor',figure_color,'Position',[800 130 70 30]);
edit_lever_threshtime_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(threshold_time_default),...
    'Position',[800 100 70 30],'Callback',@updateThreshold_manual);

% Create saving option button and save variable, nothing really necessary
% to save at the moment though so disabled
AP_labjack_lever_save = {};
save_lever_h = uicontrol(lever_figure,'Style','togglebutton',...
    'String','Save','FontSize',30, ...
    'Position',[775 15 100 70],'Callback',@saveLever,'Enable','off');


%% Callback for start button
    function start_callback(obj,event)
        
        % Check to see that timers aren't already running
        if isempty(timerfind('Running','on'));
            
            % Start timer functions
            start(getLever_timer);
            start(plotLever_timer);
            
            disp('Running...')
            
        end
    end

%% Callback for stop button
    function stop_callback(obj,event)
        
        % Check to see that timers are running
        if ~isempty(timerfind('Running','on'));
            
            % Stop and delete all timers
            stop(plotLever_timer);
            stop(getLever_timer);
            
            % This probably isn't necessary, interferes with start/stop
            %         delete(plotLever_timer);
            %         delete(getLever_timer);
            %
            %         % Close the labjack object
            %         ljudObj.Close();
            
            disp('Stopped.')
            
        end
    end


%% Function for manually updating lever trajectory threshold
    function updateThreshold_manual(obj,event)
        
        % Read values, store
        curr_restthresh = str2num(get(edit_lever_restthresh_h,'String'));
        curr_pressthresh = str2num(get(edit_lever_pressthresh_h,'String'));
        curr_threshtime = str2num(get(edit_lever_threshtime_h,'String'));
        
        lever_threshold.rest_threshold = curr_restthresh;
        lever_threshold.press_threshold = curr_pressthresh;
        lever_threshold.threshold_time = curr_threshtime;
        
        % Update the analysis buffer size
        analyze_samples = ceil(curr_threshtime)/(getLever_period*1000) + 1;
        lever_analyze_buffer = nan(analyze_samples,1);
        lever_time_analyze_buffer = nan(analyze_samples,1);

        
    end

%% Function to update the save functionality
    function saveLever(obj,event)
        [savename,savepath] = uiputfile('C:\Users\komiyama\Documents\MATLAB\ratter\SoloData\Data\experimenter\*_levertrajectory.mat');
        if savename
            save([savepath savename],'AP_labjack_lever_save');
            disp('Saved')
        else
            disp('No save destination chosen')
        end
    end

%% Timer function for reading lever
    function getLever(obj,event)
        
        ljudObj.GoOne(ljhandle);
        [ljerror, ioType, channel, curr_lever, dummyInt, dummyDbl] = ljudObj.GetFirstResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
        [ljerror, ioType, channel, curr_touch, dummyInt, dummyDbl] = ljudObj.GetNextResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
        
        % Update touch state
        touch_state = curr_touch;
        
        % Change color of strip chart depending on touch state
        switch touch_state
            case false
                set(strip_chart,'color','r')
            case true
                set(strip_chart,'color','g')
        end
        
        % Update analyzed buffers
        lever_analyze_buffer = circshift(lever_analyze_buffer,[-1 0]);
        lever_analyze_buffer(end) = curr_lever;
        
        lever_time_analyze_buffer = circshift(lever_analyze_buffer,[-1 0]);
        lever_time_analyze_buffer(end) = now;
      
        % Update plot buffers
        lever_plot_buffer = circshift(lever_plot_buffer,[-1 0]);
        lever_plot_buffer(end) = curr_lever;
        
        lever_rate = circshift(lever_rate,[-1 0]);
        lever_rate(end) = 1./getLever_period - 1./getLever_timer.InstantPeriod;
        
        % Shift output signal 
        output_signal = circshift(output_signal,[-1 0]);
        
        % Analyze lever, decide output voltage (only if touch and not 
        % within refractory period)
        
        % Define successful press: goes from rest to press threshold within
        % time limit (as defined by computer clock: can't just go on number
        % of samples because rate isn't always constant)
        successful_press = lever_analyze_buffer(end) <= lever_threshold.press_threshold && ...
            any(lever_analyze_buffer >= lever_threshold.rest_threshold) && ...
            (lever_time_analyze_buffer(find(lever_analyze_buffer(end:-1:1) >= lever_threshold.rest_threshold,1)) - ...
            lever_plot_buffer(end))*1e8 <= lever_threshold.threshold_time;
        
        if touch_state && (etime(clock,last_success_time) > refractory_period) && ...
                successful_press;
            
            % Fill output signal with ONs
            output_signal(:) = 5;
            
            % Cycle output signal to next step
            ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PUT_DAC, 0, output_signal(1), 0); 
            
            % Read out and buffer output voltage
            [ljerror, ioType, channel, curr_output, dummyInt, dummyDbl] = ljudObj.GetFirstResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
            output_buffer = circshift(output_buffer,[-1 0]);
            output_buffer(end) = curr_output;
            
            % Resubmit request for lever
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_AIN_DIFF, 5, 0, 4, 0);
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 7, 0, 0, 0);
            

            % Update last success time
            last_success_time = clock;
            
        else
            
            % Clear last value of output signal
            output_signal(end) = 0;
            
            % Cycle output signal to next step
            ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PUT_DAC, 0, output_signal(1), 0);
            
            % Read out and store output voltage
            [ljerror, ioType, channel, curr_output, dummyInt, dummyDbl] = ljudObj.GetFirstResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
            output_buffer = circshift(output_buffer,[-1 0]);
            output_buffer(end) = curr_output;
            
            % Resubmit request for lever
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_AIN_DIFF, 5, 0, 4, 0);
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 7, 0, 0, 0);
            
        end
        
    end

%% Timing function for plotting lever and sampling rate
    function plotLever(obj,event)
        set(lever_plot_h,'YData',lever_plot_buffer);
        set(output_plot_h,'YData',output_buffer);
        set(rate_plot_h,'YData',lever_rate);
        drawnow;
    end


end
