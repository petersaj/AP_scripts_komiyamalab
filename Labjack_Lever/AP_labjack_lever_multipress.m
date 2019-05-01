function AP_labjack_lever_multipress
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

trajectory_plot = subplot(4,2,[2 4]); hold on;
createTrajectory; 

output_plot = subplot(4,2,5); hold on;
xlim([0 plot_samples]);
ylim([-1 6]);
title('Output voltage')

rate_plot = subplot(4,2,7); hold on;
xlim([0 plot_samples]);
title('Sampling rate offset')

% Set definition of a press (low/high threshold, time between)
lever_trajectory.press_high_thresh = press_high_thresh;
lever_trajectory.press_low_thresh = press_low_thresh;
lever_trajectory.press_travel_time = press_travel_time;

lever_trajectory.number_presses = number_presses;
lever_trajectory.interpress_times_min = interpress_times_min;
lever_trajectory.interpress_times_max = interpress_times_max;

% Set buffer sizes for analysis, add on pre-movement samples
pre_samples = 20;
post_samples = 20;
analyze_samples = lever_trajectory.interpress_times_max/getLever_period;
buffer_samples = analyze_samples + pre_samples + post_samples;
trajectory_plot_samples = [1:buffer_samples] - pre_samples;
set(trajectory_plot,'Xlim',[trajectory_plot_samples(1) trajectory_plot_samples(end)]);

% Plot lines denoting analyzed portion
line([1 1],get(trajectory_plot,'Ylim'),'color','k','Parent',trajectory_plot);
line([analyze_samples+1 analyze_samples+1],get(trajectory_plot,'Ylim'),'color','k','Parent',trajectory_plot);

% Initialize buffers/moving flags
lever_analyze_buffer = nan(buffer_samples,1);
reward_flag = false(buffer_samples,1);
lever_plot_buffer = nan(plot_samples,1);
lever_rate = nan(plot_samples,1);
output_buffer = nan(plot_samples,1);

% Create flag for lever press (only for saving, less reliable for trigger?)
lever_press_flag = false(buffer_samples,1);

% Create a running output signal (because one sample doesn't trigger)
output_signal = zeros(5,1);

% Set initial cue state
cue_state = nan;

% Set clocks for below low / above high thresholds
last_low_threshold_time = tic;
last_press_time = tic;
last_press_intervals = nan(lever_trajectory.number_presses-1,1);

% Set last success time and refractory period (to avoid close multi-trigger)
last_success_time = tic;
refractory_period = 1; % seconds

% Create plot handles for lever and sampling rate
lever_plot_h = plot(strip_chart,lever_plot_buffer,'k');
output_plot_h = plot(output_plot,output_buffer,'k');
rate_plot_h = plot(rate_plot,lever_rate,'k');
last_successful_trajectory_h = plot(trajectory_plot, ...
    trajectory_plot_samples,lever_analyze_buffer','b','linewidth',2);

% Create start/stop buttons
start_h = uicontrol('Style','pushbutton', ...
    'Position',[550 25 100 70],'String','Start','FontSize',30, ...
    'Callback',@start_callback);

stop_h = uicontrol('Style','pushbutton', ...
    'Position',[650 25 100 70],'String','Stop','FontSize',30, ...
    'Callback',@stop_callback);

% Create saving option button and save variable, set success counter
% save variable: lever snippet, timestamp, samples with presses, intervals
success_counter = 0;
AP_labjack_lever_save = cell(0,4);
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

%% Create trajectory (in this case - presses separated by some time)
    function createTrajectory(obj,event)
        % In this incarnation: a press is defined by traveling between two
        % thresholds within a set amount of time. A number of presses with
        % defined min/max inter-press times are given, where a "press"
        % cannot be made in the middle
        
        % Set definition of a press (low/high threshold, time between)
        press_high_thresh = -1.4;
        press_low_thresh = -1.6;
        press_travel_time = 0.1; % seconds
        
        number_presses = 2;
        interpress_times_min(1) = 0.3; % seconds
        interpress_times_max(1) = 0.6; % seconds
        
        line([-100 100],[press_high_thresh press_high_thresh], ...
            'color','g','linestyle','--','Parent',trajectory_plot);
        line([-100 100],[press_low_thresh press_low_thresh], ...
            'color','r','linestyle','--','Parent',trajectory_plot);
        
        set(trajectory_plot,'YLim',[press_low_thresh - 0.5 ...
            press_high_thresh + 0.5]);

    end

%% Function for auto-updating lever trajectory

    function updateTrajectory(obj,event)
        
      % no updating at the moment
      
    end

%% Function for manually updating lever trajectory

    function updateTrajectory_manual(obj,event)     
        
        % no updating at the moment
        
    end


%% Function to save
    function saveLever(obj,event)
        
        %cd C:\Users\komiyama\Documents\MATLAB\ratter\SoloData\Data\experimenter
        [savename,savepath] = uiputfile('C:\Users\komiyama\Documents\MATLAB\ratter\SoloData\Data\experimenter\*.mat');
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
        
        % Update reward flag
        reward_flag = circshift(reward_flag,[-1 0]);
        reward_flag(end) = false;
        
        % Update lever press flag
        lever_press_flag = circshift(lever_press_flag,[-1 0]);
        lever_press_flag(end) = false;
        
        % Shift output signal 
        output_signal = circshift(output_signal,[-1 0]);
        
        % Save/plot the lever buffer if the reward flag in right position
        % (this ensures that saved/displayed buffers are aligned and the
        % same size regardless of trajectory progress)
        if reward_flag(pre_samples+1)
            % Store success info
            % save variable: lever snippet, trajectory, timestamp
            success_counter = success_counter + 1;            
            % lever snippet
            AP_labjack_lever_save{success_counter,1} = lever_analyze_buffer;
            % timestamp
            AP_labjack_lever_save{success_counter,2} = clock;
            % samples with presses
            AP_labjack_lever_save{success_counter,3} = lever_press_flag;
            % intervals
            AP_labjack_lever_save{success_counter,4} = last_press_intervals;
            
            % Plot the last successful trajectory
            set(last_successful_trajectory_h,'YData',lever_analyze_buffer);
            drawnow;
            
        end
        
        % Check if current position is below low / above high thresh
        if curr_lever < lever_trajectory.press_low_thresh
            last_low_threshold_time = tic;
        elseif curr_lever > lever_trajectory.press_high_thresh
            % if threshes crossed within time during cue, define press
            last_transthreshold_time = toc(last_low_threshold_time);
            if cue_state && ...
                    last_transthreshold_time <= lever_trajectory.press_travel_time
                last_press_intervals = circshift(last_press_intervals,[-1 0]);
                last_press_intervals(end) = toc(last_press_time);
                last_press_time = tic;
                lever_press_flag(end) = true;
            end
        end
                
        % Check if time between presses is within limit, decide output voltage 
        % (only in cue period and not within refractory period)
        if cue_state && toc(last_success_time) > refractory_period && ...
                all(last_press_intervals > lever_trajectory.interpress_times_min) && ...
                all(last_press_intervals < lever_trajectory.interpress_times_max);
            
            % Clear the press interval array
            last_press_intervals(:) = nan;
            
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
            last_success_time = tic;
            
            % Update reward flag
            reward_flag(pre_samples+post_samples-1) = true;
            
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
% This was a seperate function in hope of getting the analyzing timer
% function to run faster than this one, but it doesn't seem to matter
    function plotLever(obj,event)
        set(lever_plot_h,'YData',lever_plot_buffer);
        set(output_plot_h,'YData',output_buffer);
        set(rate_plot_h,'YData',lever_rate);
        drawnow;
    end


end
