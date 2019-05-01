%% Load, pre-process, threshold data
close all
clear all

data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\SOM_ChR2_patch';
save_filename = [data_dir filesep 'spike_latencies.mat'];
cd(data_dir);

sessions = dir(data_dir);
sessions = sessions(3:end);

if exist(save_filename,'file')
    load(save_filename)
    start_session = length(spike_latency_all);
    start_cell = length(spike_latency_all{start_session});
else
    spike_latency_all = {};
    start_session = 1;
    start_cell = 1;
end

for curr_session = start_session:length(sessions)
    
    curr_session_path = [data_dir filesep sessions(curr_session).name];
    session_dir = dir(curr_session_path);
    cell_paths = session_dir([session_dir.isdir]);
    cell_paths = cell_paths(3:end);
    
    for curr_cell = start_cell:length(cell_paths);
        
        curr_cell_path = [curr_session_path filesep cell_paths(curr_cell).name];
        cell_dir = dir(curr_cell_path);
        recording_paths = cell_dir([cell_dir.isdir]);
        recording_paths = recording_paths(3:end);
        
        for curr_recording = 1:length(recording_paths)
            
            curr_recording_path = [curr_cell_path filesep recording_paths(curr_recording).name];
            
            curr_data = AP_load_xsg_continuous(curr_recording_path);
            
            curr_drift = smooth(curr_data.channels(:,1),1000);
            curr_ephys = curr_data.channels(:,1) - curr_drift;
            curr_opto = curr_data.channels(:,2);
            
            ephys_fig = figure;
            p1 = subplot(2,1,1);
            plot(curr_data.channels(:,1),'k');
            p2 = subplot(2,1,2); hold on;
            plot(curr_ephys,'k');
            linkaxes([p1 p2],'x');
            
            title(['Session ' num2str(curr_session) ', cell ' ...
                num2str(curr_cell) ', recording ' num2str(curr_recording)]);
            
            start_sample = input('Start sample (or skip): ','s');
            switch start_sample
                case 'skip'
                    close(ephys_fig);
                    continue
                otherwise
                    start_sample = str2num(start_sample);
            end
            
            end_sample = input('End sample: ','s');
            if strcmp(end_sample,'end')
                use_samples = start_sample:size(curr_data.channels,1);
            else
                end_sample = str2num(end_sample);
                use_samples = start_sample:end_sample;
            end
            
            curr_ephys = curr_ephys(use_samples);
            curr_opto = curr_opto(use_samples);
            
            spike_thresh_validate = false;
            while ~spike_thresh_validate
                if exist('spike_plot','var')
                    delete(spike_plot);
                    clear spike_plot;
                end
                
                spike_thresh = input('Spike threshold: ');
                
                % Get spike times
                ephys_over_thresh = curr_ephys >= spike_thresh;
                ephys_under_thresh = curr_ephys < spike_thresh;
                spike_trace = double(diff(ephys_over_thresh(2:end)) == 1 & ...
                    diff(ephys_under_thresh(1:end-1)) == 0);
                spike_times = find(spike_trace);
                if isempty(spike_times)
                    disp('No spikes detected, invalid threshold')
                    continue
                end
                
                spike_plot = plot(p2,spike_times+(start_sample-1),spike_thresh,'.r');
                
                spike_check = input('Spikes ok? (y/n): ','s');
                if strcmp(spike_check,'y')
                    spike_thresh_validate = true;
                end
            end
            clear spike_plot;
            
            close(ephys_fig);
            
            % Get spike times
            ephys_over_thresh = curr_ephys >= spike_thresh;
            ephys_under_thresh = curr_ephys < spike_thresh;
            spike_trace = double(diff(ephys_over_thresh(2:end)) == 1 & ...
                diff(ephys_under_thresh(1:end-1)) == 0);
            spike_times = find(spike_trace);
            
            spike_rate = zeros(size(spike_trace));
            for curr_spike = 1:length(spike_times)-1
                spike_rate(spike_times(curr_spike):spike_times(curr_spike+1)) = ...
                    10000/(spike_times(curr_spike+1) - spike_times(curr_spike));
            end
            
            % Get opto pulse times
            opto_thresh = 2.5;
            opto_over_thresh = curr_opto >= opto_thresh;
            opto_under_thresh = curr_opto < opto_thresh;
            opto_trace = diff(opto_over_thresh(2:end)) == 1 & ...
                diff(opto_under_thresh(1:end-1)) == 0;
            opto_times = find(opto_trace);
            
            % Analyze
            
            % % Get spikes/ephys around each opto pulse
            % pre_samples = 100;
            % post_samples = 500;
            % exclude_opto_times = (opto_times - pre_samples < 1) | ...
            %     (opto_times + post_samples > length(curr_ephys));
            % opto_aligned_spikes = cell2mat(cellfun(@(x) ...
            %     spike_trace(x-pre_samples:x+post_samples)', ...
            %     num2cell(opto_times(~exclude_opto_times)),'uni',false));
            % opto_aligned_ephys = cell2mat(cellfun(@(x) ...
            %     curr_ephys(x-pre_samples:x+post_samples)', ...
            %     num2cell(opto_times(~exclude_opto_times)),'uni',false));
            
            % Get spikes for baseline + opto pulse epochs
            opto_rate = 3; %hz
            opto_rate_samples = 10000/opto_rate;
            % judge epochs by </> 2*rate
            epoch_borders = find(diff(opto_times) > opto_rate_samples*2);
            opto_onsets = unique([opto_times(1);opto_times(epoch_borders+1)]);
            opto_offsets = unique([opto_times(epoch_borders);opto_times(end)]);
            
            % get time to nearest spike for each opto pulse within epoch
            opto_epoch_pulses = mat2cell(opto_times,diff([0;epoch_borders;length(opto_times)]));
            opto_spike_distance = cellfun(@(x) arrayfun(@(y) min(spike_times(spike_times >  ...
                x(y)) - x(y))/10,1:length(x),'uni',false),opto_epoch_pulses,'uni',false);
            spike_latency = cellfun(@(x) horzcat(x{:}),opto_spike_distance,'uni',false);
            % pad spike latencies to all be the same size
            max_pulses = max(cellfun(@length,spike_latency));
            spike_latency_pad = cell2mat(cellfun(@(x) padarray(x, ...
                [0 max_pulses-length(x)],nan,'post'),spike_latency,'uni',false));
            
            
            % pre_opto = nanmedian(opto_offsets - opto_onsets);
            % post_opto = 30;
            %
            % % sometimes different number of pulses, pad by max time
            % max_epoch_time = max(opto_offsets(~exclude_epochs) - opto_onsets(~exclude_epochs)) + pre_opto + post_opto + 1;
            %
            % % don't use epochs that don't have a baseline or baseline's worth post-time
            % exclude_epochs = (opto_onsets - pre_opto < 1) | ...
            %     (opto_offsets + pre_opto > length(curr_ephys));
            %
            % opto_epoch_aligned_spikes = cell2mat(cellfun(@(x,y) ...
            %     padarray(spike_trace(x-pre_opto:y+post_opto)',[0 max_epoch_time - ...
            %     length(spike_trace(x-pre_opto:y+post_opto))],nan,'post'),num2cell(opto_onsets(~exclude_epochs)), ...
            %     num2cell(opto_offsets(~exclude_epochs)),'uni',false));
            % opto_epoch_aligned_spikerate = cell2mat(cellfun(@(x,y) ...
            %     padarray(spike_rate(x-pre_opto:y+post_opto)',[0 max_epoch_time - ...
            %     length(spike_rate(x-pre_opto:y+post_opto))],nan,'post'),num2cell(opto_onsets(~exclude_epochs)), ...
            %     num2cell(opto_offsets(~exclude_epochs)),'uni',false));
            % opto_epoch_aligned_ephys = cell2mat(cellfun(@(x,y) ...
            %     padarray(curr_ephys(x-pre_opto:y+post_opto)',[0 max_epoch_time - ...
            %     length(curr_ephys(x-pre_opto:y+post_opto))],nan,'post'),num2cell(opto_onsets(~exclude_epochs)), ...
            %     num2cell(opto_offsets(~exclude_epochs)),'uni',false));
            
            % Save spike latencies (in ms)
            spike_latency_all{curr_session}{curr_cell}{curr_recording} = horzcat(spike_latency{:});
            
            save(save_filename,'spike_latency_all');
        end
    end
end

%% To re-do a given cell/session
close all
clear all

curr_session = 5;
curr_cell = 1;

data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\SOM_ChR2_patch';
save_filename = [data_dir filesep 'spike_latencies.mat'];
load(save_filename)

% Define session and cells
sessions = dir(data_dir);
sessions = sessions(3:end);

curr_session_path = [data_dir filesep sessions(curr_session).name];
session_dir = dir(curr_session_path);
cell_paths = session_dir([session_dir.isdir]);
cell_paths = cell_paths(3:end);

% Get recordings, loop through
curr_cell_path = [curr_session_path filesep cell_paths(curr_cell).name];
cell_dir = dir(curr_cell_path);
recording_paths = cell_dir([cell_dir.isdir]);
recording_paths = recording_paths(3:end);

for curr_recording = 1:length(recording_paths)
    
    curr_recording_path = [curr_cell_path filesep recording_paths(curr_recording).name];
    
    curr_data = AP_load_xsg_continuous(curr_recording_path);
    
    curr_drift = smooth(curr_data.channels(:,1),1000);
    curr_ephys = curr_data.channels(:,1) - curr_drift;
    curr_opto = curr_data.channels(:,2);
    
    ephys_fig = figure;
    p1 = subplot(2,1,1);
    plot(curr_data.channels(:,1),'k');
    p2 = subplot(2,1,2); hold on;
    plot(curr_ephys,'k');
    linkaxes([p1 p2],'x');
    
    title(['Session ' num2str(curr_session) ', cell ' ...
        num2str(curr_cell) ', recording ' num2str(curr_recording)]);
    
    start_sample = input('Start sample (or skip): ','s');
    switch start_sample
        case 'skip'
            spike_latency_all{curr_session}{curr_cell}{curr_recording} = [];
            close(ephys_fig);
            continue
        otherwise
            start_sample = str2num(start_sample);
    end
    
    end_sample = input('End sample: ','s');
    if strcmp(end_sample,'end')
        use_samples = start_sample:size(curr_data.channels,1);
    else
        end_sample = str2num(end_sample);
        use_samples = start_sample:end_sample;
    end
    
    curr_ephys = curr_ephys(use_samples);
    curr_opto = curr_opto(use_samples);
    
    spike_thresh_validate = false;
    while ~spike_thresh_validate
        if exist('spike_plot','var')
            delete(spike_plot);
            clear spike_plot;
        end
        
        spike_thresh = input('Spike threshold: ');
        
        % Get spike times
        ephys_over_thresh = curr_ephys >= spike_thresh;
        ephys_under_thresh = curr_ephys < spike_thresh;
        spike_trace = double(diff(ephys_over_thresh(2:end)) == 1 & ...
            diff(ephys_under_thresh(1:end-1)) == 0);
        spike_times = find(spike_trace);
        if isempty(spike_times)
            disp('No spikes detected, invalid threshold')
            continue
        end
        
        spike_plot = plot(p2,spike_times+(start_sample-1),spike_thresh,'.r');
        
        spike_check = input('Spikes ok? (y/n): ','s');
        if strcmp(spike_check,'y')
            spike_thresh_validate = true;
        end
    end
    clear spike_plot;
    
    close(ephys_fig);
    
    % Get spike times
    ephys_over_thresh = curr_ephys >= spike_thresh;
    ephys_under_thresh = curr_ephys < spike_thresh;
    spike_trace = double(diff(ephys_over_thresh(2:end)) == 1 & ...
        diff(ephys_under_thresh(1:end-1)) == 0);
    spike_times = find(spike_trace);
    
    spike_rate = zeros(size(spike_trace));
    for curr_spike = 1:length(spike_times)-1
        spike_rate(spike_times(curr_spike):spike_times(curr_spike+1)) = ...
            10000/(spike_times(curr_spike+1) - spike_times(curr_spike));
    end
    
    % Get opto pulse times
    opto_thresh = 2.5;
    opto_over_thresh = curr_opto >= opto_thresh;
    opto_under_thresh = curr_opto < opto_thresh;
    opto_trace = diff(opto_over_thresh(2:end)) == 1 & ...
        diff(opto_under_thresh(1:end-1)) == 0;
    opto_times = find(opto_trace);
    
    % Get spikes for baseline + opto pulse epochs
    opto_rate = 3; %hz
    opto_rate_samples = 10000/opto_rate;
    % judge epochs by </> 2*rate
    epoch_borders = find(diff(opto_times) > opto_rate_samples*2);
    opto_onsets = unique([opto_times(1);opto_times(epoch_borders+1)]);
    opto_offsets = unique([opto_times(epoch_borders);opto_times(end)]);
    
    % get time to nearest spike for each opto pulse within epoch
    opto_epoch_pulses = mat2cell(opto_times,diff([0;epoch_borders;length(opto_times)]));
    opto_spike_distance = cellfun(@(x) arrayfun(@(y) min(spike_times(spike_times >  ...
        x(y)) - x(y))/10,1:length(x),'uni',false),opto_epoch_pulses,'uni',false);
    spike_latency = cellfun(@(x) horzcat(x{:}),opto_spike_distance,'uni',false);
    
    % Save spike latencies (in ms)
    spike_latency_all{curr_session}{curr_cell}{curr_recording} = horzcat(spike_latency{:});
    
end
save(save_filename,'spike_latency_all');


%% Plot latency histogram after manual thresholding / checking

clear all
data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\SOM_ChR2_patch';
save_filename = [data_dir filesep 'spike_latencies.mat'];
load(save_filename)

% get rid of empty cells (what happened there?)
spike_latency_use = cellfun(@(x) x(~cellfun(@isempty,x)), spike_latency_all,'uni',false);

% this is just to manually check that the recordings in cell are similar
epoch_median_latency = cellfun(@(x) cellfun(@(x) cellfun(@nanmedian,x),x, ...
    'uni',false),spike_latency_use,'uni',false);

% concatenate epochs, get median for cell
spike_latency_allrec = cellfun(@(x) cellfun(@(x) ...
    horzcat(x{:}),x,'uni',false),spike_latency_use,'uni',false);

spike_latency_allcell = horzcat(spike_latency_allrec{:});

% get rid of spikes > 100 ms (definitely not opto-elicited)
spike_latency_optospike = cellfun(@(x) x(x < 100),spike_latency_allcell,'uni',false);

spike_latency_median = cellfun(@nanmedian,spike_latency_optospike);
spike_latency_std = cellfun(@nanstd,spike_latency_optospike);

figure;
errorbar(spike_latency_median,spike_latency_std,'.k','linewidth',2);


%% Make figure of the 10 minute cell / 3 Hz

close all
clear all

curr_session = 3;
curr_cell = 2;
curr_recording = 5;

data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\SOM_ChR2_patch';
save_filename = [data_dir filesep 'spike_latencies.mat'];
load(save_filename)

% Define session and cells
sessions = dir(data_dir);
sessions = sessions(3:end);

curr_session_path = [data_dir filesep sessions(curr_session).name];
session_dir = dir(curr_session_path);
cell_paths = session_dir([session_dir.isdir]);
cell_paths = cell_paths(3:end);

% Get recordings, loop through
curr_cell_path = [curr_session_path filesep cell_paths(curr_cell).name];
cell_dir = dir(curr_cell_path);
recording_paths = cell_dir([cell_dir.isdir]);
recording_paths = recording_paths(3:end);

curr_recording_path = [curr_cell_path filesep recording_paths(curr_recording).name];

curr_data = AP_load_xsg_continuous(curr_recording_path);

curr_drift = smooth(curr_data.channels(:,1),1000);
curr_ephys = curr_data.channels(:,1) - curr_drift;
curr_opto = curr_data.channels(:,2);

% Get spike times
spike_thresh = 1;
ephys_over_thresh = curr_ephys >= spike_thresh;
ephys_under_thresh = curr_ephys < spike_thresh;
spike_trace = double(diff(ephys_over_thresh(2:end)) == 1 & ...
    diff(ephys_under_thresh(1:end-1)) == 0);
spike_times = find(spike_trace);

spike_rate = zeros(size(spike_trace));
for curr_spike = 1:length(spike_times)-1
    spike_rate(spike_times(curr_spike):spike_times(curr_spike+1)) = ...
        10000/(spike_times(curr_spike+1) - spike_times(curr_spike));
end

% Get opto pulse times
opto_thresh = 2.5;
opto_over_thresh = curr_opto >= opto_thresh;
opto_under_thresh = curr_opto < opto_thresh;
opto_trace = diff(opto_over_thresh(2:end)) == 1 & ...
    diff(opto_under_thresh(1:end-1)) == 0;
opto_times = find(opto_trace);

% Get spikes for baseline + opto pulse epochs
opto_rate = 3; %hz
opto_rate_samples = 10000/opto_rate;
% judge epochs by </> 2*rate
epoch_borders = find(diff(opto_times) > opto_rate_samples*2);
opto_onsets = unique([opto_times(1);opto_times(epoch_borders+1)]);
opto_offsets = unique([opto_times(epoch_borders);opto_times(end)]);

% get time to nearest spike for each opto pulse within epoch
opto_epoch_pulses = mat2cell(opto_times,diff([0;epoch_borders;length(opto_times)]));
spike_latency = cell2mat(cellfun(@(x) arrayfun(@(y) min(spike_times(spike_times >  ...
    x(y)) - x(y))/10,1:length(x)),opto_epoch_pulses,'uni',false));


figure;

t_start = 50000;
t1 = t_start:t_start+2*10000;
subplot(1,2,1); hold on;
plot(curr_opto(t1) > 2.5,'k')
xlim([0 t1(end)-t1(1)])
ylim([-0.5 1.5])

t1_spikes = spike_times(spike_times >= t1(1) & spike_times <= t1(end))-t1(1)+1;
for i = 1:length(t1_spikes);
    line(repmat(t1_spikes(i),2,1),[0.5 1],'color','r','linewidth',2);
end

% 10 mins later (and a little to center it)
fudge_factor = 0;
t2 = t1(end)+10*60*10000 + fudge_factor : t1(end)+10*60*10000 + 2*10000 + fudge_factor;
subplot(1,2,2);
plot(curr_opto(t2) > 2.5,'k')
xlim([0 t2(end)-t2(1)])
ylim([-0.5 1.5])

t2_spikes = spike_times(spike_times >= t2(1) & spike_times <= t2(end))-t2(1)+1;
for i = 1:length(t2_spikes);
    line(repmat(t2_spikes(i),2,1),[0.5 1],'color','r','linewidth',2);
end


%% Make figure of multiple trials with raw traces

close all
clear all

curr_session = 4;
curr_cell = 6;
curr_recording = 1;

data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\SOM_ChR2_patch';
save_filename = [data_dir filesep 'spike_latencies.mat'];
load(save_filename)

% Define session and cells
sessions = dir(data_dir);
sessions = sessions(3:end);

curr_session_path = [data_dir filesep sessions(curr_session).name];
session_dir = dir(curr_session_path);
cell_paths = session_dir([session_dir.isdir]);
cell_paths = cell_paths(3:end);

% Get recordings, loop through
curr_cell_path = [curr_session_path filesep cell_paths(curr_cell).name];
cell_dir = dir(curr_cell_path);
recording_paths = cell_dir([cell_dir.isdir]);
recording_paths = recording_paths(3:end);

curr_recording_path = [curr_cell_path filesep recording_paths(curr_recording).name];

curr_data = AP_load_xsg_continuous(curr_recording_path);

curr_drift = smooth(curr_data.channels(:,1),1000);
curr_ephys = curr_data.channels(:,1) - curr_drift;
curr_opto = curr_data.channels(:,2);

% Get opto pulse times
opto_thresh = 2.5;
opto_over_thresh = curr_opto >= opto_thresh;
opto_under_thresh = curr_opto < opto_thresh;
opto_trace = diff(opto_over_thresh(2:end)) == 1 & ...
    diff(opto_under_thresh(1:end-1)) == 0;
opto_times = find(opto_trace);

% Get spikes for baseline + opto pulse epochs
opto_rate = 3; %hz
opto_rate_samples = 10000/opto_rate;
% judge epochs by </> 2*rate
epoch_borders = find(diff(opto_times) > opto_rate_samples*2);
opto_onsets = unique([opto_times(1);opto_times(epoch_borders+1)]);
opto_offsets = unique([opto_times(epoch_borders);opto_times(end)]);


pre_opto = nanmedian(opto_offsets - opto_onsets);
post_opto = 0;

% don't use epochs that don't have a baseline or baseline's worth post-time
exclude_epochs = (opto_onsets - pre_opto < 1) | ...
    (opto_offsets + pre_opto > length(curr_ephys));

% sometimes different number of pulses, pad by max time
max_epoch_time = max(opto_offsets(~exclude_epochs) - opto_onsets(~exclude_epochs)) + pre_opto + post_opto + 1;

opto_epoch_aligned_ephys = cell2mat(cellfun(@(x,y) ...
    padarray(curr_ephys(x-pre_opto:y+post_opto)',[0 max_epoch_time - ...
    length(curr_ephys(x-pre_opto:y+post_opto))],nan,'post'),num2cell(opto_onsets(~exclude_epochs)), ...
    num2cell(opto_offsets(~exclude_epochs)),'uni',false));

% Plot traces for a few epochs
figure;
p1 = subplot(10,1,1);
plot_opto_epoch = opto_onsets(find(~exclude_epochs,1)) - pre_opto: ...
    opto_offsets(find(~exclude_epochs,1)) + post_opto;
plot(curr_opto(plot_opto_epoch) > 2.5,'k');
axis off

p2 = subplot(10,1,2:10); hold on;
for i = 1:10;
    plot(opto_epoch_aligned_ephys(i,:) + 500*i,'k');
end
linkaxes([p1 p2],'x');

% Plot one pulse for all trials
figure;
plot_time = pre_opto-100:pre_opto + 300;
p1 = subplot(10,1,1);
plot(curr_opto(plot_opto_epoch(plot_time)) > 2.5,'k');
axis off

p2 = subplot(10,1,2:10); hold on;
for i = 1:size(opto_epoch_aligned_ephys,1);
    plot(opto_epoch_aligned_ephys(i,plot_time) + 500*i,'k');
end
linkaxes([p1 p2],'x');








