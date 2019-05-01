function AP_labjack_lever_threshold

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

% Get, group current trajectory information
lever_trajectory.threshold = y_threshold;
lever_trajectory.sign = y_sign;
% Set default starting samples to analyze (aka progress)
lever_trajectory.progress = 12;
% Set default number of hits to progress
lever_trajectory.level_up = 5;
lever_trajectory.level_progress = 1;

% Set buffer sizes for analysis, add on pre-movement samples
pre_samples = 10;
analyze_samples = length(lever_trajectory.threshold) + pre_samples;
set(trajectory_plot,'Xlim',[ -pre_samples+1 length(lever_trajectory.threshold)]);

% Plot line at zero point where analzyed movement begins
line([1 1],get(trajectory_plot,'Ylim'),'color','k','Parent',trajectory_plot);

% Initialize buffers/moving flags
lever_analyze_buffer = nan(analyze_samples,1);
reward_flag = false(analyze_samples,1);
lever_plot_buffer = nan(plot_samples,1);
lever_rate = nan(plot_samples,1);
output_buffer = nan(plot_samples,1);
% Create a running output signal (because one sample doesn't trigger)
output_signal = zeros(5,1);

% Set initial cue state
cue_state = NaN;

% Set last success time and refractory period (to avoid close multi-trigger)
last_success_time = clock;
refractory_period = 1; % seconds

% Create plot handles for lever and sampling rate
lever_plot_h = plot(strip_chart,lever_plot_buffer,'k');
output_plot_h = plot(output_plot,output_buffer,'k');
rate_plot_h = plot(rate_plot,lever_rate,'k');
last_successful_trajectory_h = plot(trajectory_plot, ...
    -pre_samples+1:length(lever_trajectory.threshold), ...
    lever_analyze_buffer','b','linewidth',2);

% Create start/stop buttons
start_h = uicontrol('Style','pushbutton', ...
    'Position',[550 25 100 70],'String','Start','FontSize',30, ...
    'Callback',@start_callback);

stop_h = uicontrol('Style','pushbutton', ...
    'Position',[650 25 100 70],'String','Stop','FontSize',30, ...
    'Callback',@stop_callback);

% Create saving option button and save variable, set success counter
% save variable: lever snippet, trajectory, sign, progress, timestamp
success_counter = 0;
AP_labjack_lever_save = cell(0,5);
save_lever_h = uicontrol(lever_figure,'Style','togglebutton',...
    'String','Save','FontSize',30, ...
    'Position',[750 25 100 70],'Callback',@saveLever);

% Create editable fields for updating trajectory with success
figure_color = get(lever_figure,'Color');

uicontrol(lever_figure,'Style','text','String','Progress', ...
    'BackgroundColor',figure_color,'Position',[545 130 70 30]);
edit_lever_trajectory_progress_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.progress),...
    'Position',[545 110 70 30],'Callback',@updateTrajectory_manual);

uicontrol(lever_figure,'Style','text','String','Level up', ...
    'BackgroundColor',figure_color,'Position',[625 130 70 30]);
edit_lever_trajectory_level_up_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.level_up),...
    'Position',[625 110 70 30],'Callback',@updateTrajectory_manual);

uicontrol(lever_figure,'Style','text','String','Level progress', ...
    'BackgroundColor',figure_color,'Position',[700 130 80 30]);
edit_lever_trajectory_level_progress_h = uicontrol(lever_figure,'Style','edit',...
    'String',num2str(lever_trajectory.level_progress),...
    'Position',[705 110 70 30],'Callback',@updateTrajectory_manual);

% Draw line for current trajectory progress
trajectory_progress_h = line(repmat(lever_trajectory.progress,1,2), ...
    get(trajectory_plot,'Ylim'),'Parent', ...
    trajectory_plot,'color','k','linestyle','--','linewidth',2);

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

%% Create trajectory (in this case - series of thresholds)
    function createTrajectory(obj,event)
        
        thresh_sign{1} = -1; % (+ for over/push, - for under/pull)
        thresh{1} = -1.7; % in V
        thresh_time{1} = 0.1; % in sec
        thresh_transition_time{1} = 0.1; % in sec
        thresh_samples{1} = round(thresh_time{1}/getLever_period);
        thresh_transition_samples{1} = round(thresh_transition_time{1}/getLever_period);
        thresh_trajectory{1} = [thresh{1}*ones(1,thresh_samples{1}) ...
            -Inf*thresh_sign{1}*ones(1,thresh_transition_samples{1})];
        thresh_trajectory_sign{1} = thresh_sign{1}* ...
            ones(1,thresh_samples{1}+thresh_transition_samples{1});
        
        thresh_sign{2} = 1; % (+ for over/push, - for under/pull)
        thresh{2} = -1.4; % in V
        thresh_time{2} = 0.1; % in sec
        thresh_transition_time{2} = 0.1; % in sec
        thresh_samples{2} = round(thresh_time{2}/getLever_period);
        thresh_transition_samples{2} = round(thresh_transition_time{2}/getLever_period);
        thresh_trajectory{2} = [thresh{2}*ones(1,thresh_samples{2}) ...
            -Inf*thresh_sign{2}*ones(1,thresh_transition_samples{2})];
        thresh_trajectory_sign{2} = thresh_sign{2}* ...
            ones(1,thresh_samples{2}+thresh_transition_samples{2});
        
        thresh_sign{3} = -1; % (+ for over/push, - for under/pull)
        thresh{3} = -1.7; % in V
        thresh_time{3} = 0.1; % in sec
        thresh_transition_time{3} = 0.1; % in sec
        thresh_samples{3} = round(thresh_time{3}/getLever_period);
        thresh_transition_samples{3} = round(thresh_transition_time{3}/getLever_period);
        thresh_trajectory{3} = [thresh{3}*ones(1,thresh_samples{3}) ...
            -Inf*thresh_sign{3}*ones(1,thresh_transition_samples{3})];
        thresh_trajectory_sign{3} = thresh_sign{3}* ...
            ones(1,thresh_samples{3}+thresh_transition_samples{3});
        
        thresh_sign{4} = 1; % (+ for over/push, - for under/pull)
        thresh{4} = -1.4; % in V
        thresh_time{4} = 0.1; % in sec
        thresh_transition_time{4} = 0; % in sec
        thresh_samples{4} = round(thresh_time{4}/getLever_period);
        thresh_transition_samples{4} = round(thresh_transition_time{4}/getLever_period);
        thresh_trajectory{4} = [thresh{4}*ones(1,thresh_samples{4}) ...
            -Inf*thresh_sign{4}*ones(1,thresh_transition_samples{4})];
        thresh_trajectory_sign{4} = thresh_sign{4}* ...
            ones(1,thresh_samples{4}+thresh_transition_samples{4});
        
        y_threshold = horzcat(thresh_trajectory{:});
        y_sign = horzcat(thresh_trajectory_sign{:});
     
        plot(trajectory_plot,find(y_sign == 1),y_threshold(y_sign == 1),'g','linewidth',2); 
        plot(trajectory_plot,find(y_sign == -1),y_threshold(y_sign == -1),'r','linewidth',2); 
        
        lever_threshold_h = plot(trajectory_plot,y_threshold','--k');
        ylim([min(y_threshold(~isinf(y_threshold))) - 0.8; max(y_threshold(~isinf(y_threshold))) + 0.8])
        xlim([0 length(y_threshold)]);
        title('Thresholds')
    end

%% Function for auto-updating lever trajectory

    function updateTrajectory(obj,event)
        
        % Increment trajectory progress by 1 sample if level progressed
        if lever_trajectory.level_progress == lever_trajectory.level_up
            
            lever_trajectory.level_progress = 1;
            lever_trajectory.progress = lever_trajectory.progress+1;
            set(trajectory_progress_h,'XData',repmat(lever_trajectory.progress,1,2));
            drawnow;
            
            % Update text box
            set(edit_lever_trajectory_progress_h,'String',num2str(lever_trajectory.progress));
            set(edit_lever_trajectory_level_progress_h,'String',num2str(lever_trajectory.level_progress));
            
        else
            
            % Update level progress
            lever_trajectory.level_progress = lever_trajectory.level_progress + 1;
            % Update text box
            set(edit_lever_trajectory_level_progress_h,'String',num2str(lever_trajectory.level_progress));
            
        end
    end

%% Function for manually updating lever trajectory

    function updateTrajectory_manual(obj,event)     
        
        curr_progress = str2num(get(edit_lever_trajectory_progress_h,'String'));
        curr_level_up = str2num(get(edit_lever_trajectory_level_up_h,'String'));
        curr_level_progress = str2num(get(edit_lever_trajectory_level_progress_h,'String'));
        
        if curr_progress > length(lever_trajectory.threshold);
            curr_progress = length(lever_trajectory.threshold);
        elseif curr_progress < 2
            curr_progress = 2;
        end
        
        lever_trajectory.progress = curr_progress;
        lever_trajectory.level_up = curr_level_up;
        lever_trajectory.level_progress = curr_level_progress;
        
        set(trajectory_progress_h,'XData',repmat(lever_trajectory.progress,1,2));
        drawnow;       
        
        % Update text box
        set(edit_lever_trajectory_progress_h,'String',num2str(lever_trajectory.progress));
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
            % threshold
            AP_labjack_lever_save{success_counter,2} = lever_trajectory.threshold;
            % sign
            AP_labjack_lever_save{success_counter,3} = lever_trajectory.sign;
            % progress
            AP_labjack_lever_save{success_counter,4} = lever_trajectory.progress;
            % timestamp
            AP_labjack_lever_save{success_counter,5} = clock;
            
            % Plot the last successful trajectory
            set(last_successful_trajectory_h,'YData',lever_analyze_buffer);
            drawnow;
            
            % Progress the lever trajectory if applicable
            if lever_trajectory.progress < length(lever_trajectory.threshold)
                updateTrajectory;
            end
        end
        
        % Analyze lever amount dictated by progress, decide output voltage 
        % (only in cue period and not within refractory period)
        if cue_state && (etime(clock,last_success_time) > refractory_period) && ...
                all(lever_analyze_buffer(end-lever_trajectory.progress+1:end)'.* ...
                lever_trajectory.sign(1:lever_trajectory.progress) > ...
                lever_trajectory.threshold(1:lever_trajectory.progress).* ...
                lever_trajectory.sign(1:lever_trajectory.progress));
                       
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
            
            % Update reward flag
            reward_flag(end-lever_trajectory.progress+1) = true;
            
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
