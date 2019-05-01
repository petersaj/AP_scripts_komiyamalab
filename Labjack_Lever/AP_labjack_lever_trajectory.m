function AP_labjack_lever_trajectory

clear all
clc

%% Initialize LabJack
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
% FIO7 should be set to digital by default (cue state)

% Execute all queued commands
ljudObj.GoOne(ljhandle);

%Request a reading from FIO5-FIO4 (lever)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_AIN_DIFF, 5, 0, 4, 0);

% Request a reading from FIO7 (cue state)
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

% Start acquisition/plotting timers
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

% cue_state_plot = subplot(4,2,3); hold on;
% xlim([0 plot_samples]);
% ylim([-1 6]);
% title('Cue state');

trajectory_plot = subplot(4,2,[2 4]); hold on;
createTrajectory; 

output_plot = subplot(4,2,5); hold on;
xlim([0 plot_samples]);
ylim([-1 6]);
title('Output voltage')

rate_plot = subplot(4,2,7); hold on;
xlim([0 plot_samples]);
title('Sampling rate offset')

% Set buffer sizes for analysis, initialize
analyze_samples = length(y_trajectory);

% Get, group current trajectory information
lever_trajectory.trajectory = y_trajectory;
lever_trajectory.threshold = y_threshold;
lever_trajectory.lower_limit = y_lower_limit;
lever_trajectory.leeway = y_leeway;
lever_trajectory.plot_h = lever_threshold_h;
% Set default auto-narrow amount for threshold
lever_trajectory.auto_narrow = 0.01;
% Set default minimum leeway
lever_trajectory.minleeway = 0.1;

lever_analyze_buffer = nan(analyze_samples,1);
lever_plot_buffer = nan(plot_samples,1);
lever_rate = nan(plot_samples,1);
output_buffer = nan(plot_samples,1);

% Set initial cue state
cue_state = NaN;

% Set last success time and refractory period (in seconds)
last_success_time = clock;
refractory_period = 1;

% Create plot handles for lever and sampling rate
lever_plot_h = plot(strip_chart,lever_plot_buffer,'k');
output_plot_h = plot(output_plot,output_buffer,'k');
rate_plot_h = plot(rate_plot,lever_rate,'k');
last_successful_trajectory_h = plot(trajectory_plot, ...
    lever_analyze_buffer','b','linewidth',2);

% Create start/stop buttons
start_h = uicontrol('Style','pushbutton', ...
    'Position',[550 25 100 70],'String','Start','FontSize',30, ...
    'Callback',@start_callback);

stop_h = uicontrol('Style','pushbutton', ...
    'Position',[650 25 100 70],'String','Stop','FontSize',30, ...
    'Callback',@stop_callback);

% Create editable fields for lever threshold/floor/auto-narrow
figure_color = get(lever_figure,'Color');

uicontrol(lever_figure,'Style','text','String','Leeway', ...
    'BackgroundColor',figure_color,'Position',[545 130 70 30]);
edit_lever_threshold_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.leeway),...
    'Position',[545 110 70 30],'Callback',@updateTrajectory_manual);

uicontrol(lever_figure,'Style','text','String','Min leeway', ...
    'BackgroundColor',figure_color,'Position',[625 130 70 30]);
edit_lever_minthresh_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.minleeway),...
    'Position',[625 110 70 30],'Callback',@updateTrajectory_manual);

uicontrol(lever_figure,'Style','text','String','Floor', ...
    'BackgroundColor',figure_color,'Position',[705 130 70 30]);
edit_lever_floor_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.lower_limit),...
    'Position',[705 110 70 30],'Callback',@updateTrajectory_manual);

uicontrol(lever_figure,'Style','text','String','Auto-narrow', ...
    'BackgroundColor',figure_color,'Position',[785 130 70 30]);
edit_lever_autonarrow_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.auto_narrow),...
    'Position',[785 110 70 30],'Callback',@updateTrajectory_manual);

% Create saving option button and save variable, set success counter
% save variable: lever snippet, trajectory, timestamp
success_counter = 0;
AP_labjack_lever_save = cell(0,3);
save_lever_h = uicontrol(lever_figure,'Style','togglebutton',...
    'String','Save','FontSize',30, ...
    'Position',[750 25 100 70],'Callback',@saveLever);


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

%% Create trajectory
    function createTrajectory(obj,event)
        
        t = 0:getLever_period:2;
        y1_amp = 0.4;
        y2_amp = 0.8;
        y1 = 0.4*-sin(t*2*pi)-y1_amp;
        y2 = 0.8*sin(t*4*pi)-y2_amp;
        y_combine = y1+y2;
        
        x = t-median(t);
        mean = 0;
        sigma = 0.3;
        y_damper = (1./sigma*sqrt(2*pi))*exp(-(x-mean).^2/(2*sigma^2));
        y_damper_norm = y_damper./max(y_damper);
        
        y_gain = -0.2;
        
        y_offset = -1.8;
        
        y_trajectory = y_combine.*y_damper_norm*y_gain+y_offset;
        plot(trajectory_plot,y_trajectory,'k','linewidth',2);
        
        % define thresholds for reward
        y_lower_limit = y_offset - 0.2;
        y_leeway = 0.3;
        y_threshold = [y_trajectory-y_leeway;y_trajectory+y_leeway];
        y_threshold(y_threshold < y_lower_limit) = y_lower_limit;
        hold on
        lever_threshold_h = plot(trajectory_plot,y_threshold','--k');
        ylim([min(y_threshold(:)) - 0.2; max(y_threshold(:)) + 0.2])
        title('Trajectory')
    end


%% Function for auto-updating lever trajectory threshold

    function updateTrajectory(obj,event)
        lever_trajectory.leeway = lever_trajectory.leeway - ...
            lever_trajectory.auto_narrow;
        
        lever_trajectory.trajectory = y_trajectory;
        lever_trajectory.threshold = ...
            [lever_trajectory.trajectory - lever_trajectory.leeway; ...
            lever_trajectory.trajectory + lever_trajectory.leeway];
        lever_trajectory.threshold(lever_trajectory.threshold < ...
            lever_trajectory.lower_limit) = lever_trajectory.lower_limit;
        
        set(lever_trajectory.plot_h(1),'YData',lever_trajectory.threshold(1,:));
        set(lever_trajectory.plot_h(2),'YData',lever_trajectory.threshold(2,:));
        drawnow;       
        
        % update threshold text box
        set(edit_lever_threshold_h,'String',num2str(lever_trajectory.leeway));
        
    end

%% Function for manually updating lever trajectory threshold
    function updateTrajectory_manual(obj,event)
        
        curr_threshold = str2num(get(edit_lever_threshold_h,'String'));
        curr_floor = str2num(get(edit_lever_floor_h,'String'));
        curr_autonarrow = str2num(get(edit_lever_autonarrow_h,'String'));
        
        lever_trajectory.lower_limit = curr_floor;
        lever_trajectory.leeway = curr_threshold;
        lever_trajectory.auto_narrow = curr_autonarrow;
        
        lever_trajectory.trajectory = y_trajectory;
        lever_trajectory.threshold = ...
            [lever_trajectory.trajectory - lever_trajectory.leeway; ...
            lever_trajectory.trajectory + lever_trajectory.leeway];
        lever_trajectory.threshold(lever_trajectory.threshold < ...
            lever_trajectory.lower_limit) = lever_trajectory.lower_limit;
        
        set(lever_trajectory.plot_h(1),'YData',lever_trajectory.threshold(1,:));
        set(lever_trajectory.plot_h(2),'YData',lever_trajectory.threshold(2,:));
        drawnow;       
    end

%% Function to save
    function saveLever(obj,event)
        %cd C:\Users\komiyama\Documents\MATLAB\ratter\SoloData\Data\experimenter
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
        [ljerror, ioType, channel, curr_cue, dummyInt, dummyDbl] = ljudObj.GetNextResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
        
        % Update cue state
        cue_state = curr_cue;
        
        % Change color of strip chart depending on cue state
        switch cue_state
            case false
                set(strip_chart,'color','r')
            case true
                set(strip_chart,'color','g')
        end
        
        % Update analyzed buffer
        lever_analyze_buffer = circshift(lever_analyze_buffer,[-1 0]);
        lever_analyze_buffer(end) = curr_lever;
        
        lever_rate = circshift(lever_rate,[-1 0]);
        lever_rate(end) = 1./getLever_period - 1./getLever_timer.InstantPeriod;
        
        % Update plot buffer
        lever_plot_buffer = circshift(lever_plot_buffer,[-1 0]);
        lever_plot_buffer(end) = curr_lever;
        
        % Analyze lever, decide output voltage (only in cue period and not 
        % within refractory period)
        if cue_state && (etime(clock,last_success_time) > refractory_period) && ...
                all(lever_analyze_buffer > lever_trajectory.threshold(1,:)' & ...
                lever_analyze_buffer < lever_trajectory.threshold(2,:)');
            
            % Turn ON output voltage
            ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PUT_DAC, 0, 5, 0);
            
            % Read out and buffer output voltage
            [ljerror, ioType, channel, curr_output, dummyInt, dummyDbl] = ljudObj.GetFirstResult(ljhandle, ioDummy, chanDummy, 0, 0, 0);
            output_buffer = circshift(output_buffer,[-1 0]);
            output_buffer(end) = curr_output;
            
            % Store success info
            % save variable: lever snippet, trajectory, timestamp
            success_counter = success_counter + 1;
            % lever snippet
            AP_labjack_lever_save{success_counter,1} = lever_analyze_buffer;
            % threshold
            AP_labjack_lever_save{success_counter,2} = lever_trajectory.threshold;
            % timestamp
            AP_labjack_lever_save{success_counter,3} = clock;
            
            % Resubmit request for lever
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_AIN_DIFF, 5, 0, 4, 0);
            ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.GET_DIGITAL_BIT, 7, 0, 0, 0);
            
            % Plot the last successful trajectory
            set(last_successful_trajectory_h,'YData',lever_analyze_buffer);
            drawnow;
            
            % Auto-update trajectory (if minimum hasn't been reached)
            if lever_trajectory.leeway > lever_trajectory.minleeway
                updateTrajectory;
            end
            
            % Update last success time
            last_success_time = clock;
            
        else
            
            % Turn OFF output voltage
            ljudObj.ePut(ljhandle, LabJack.LabJackUD.IO.PUT_DAC, 0, 0, 0);
            
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
