figure;hold on;
for i = 1:length(data.coincidence_indx_tk_all);
    plot(data.coincidence_indx_tk_all{i},'.','color',[i/length(data.coincidence_indx_tk_all) 0 0]);
end

%%
num_rep = 1000;
ci = zeros(2,length(data.coincidence_indx_tk_all));
for j = 1:size(data.coincidence_indx_tk_all,2)
    ci(:,j) = bootci(num_rep,@nanmean,data.coincidence_indx_tk_all{j});
    disp(['Bootstrapping ' num2str(j/length(data.coincidence_indx_tk_all))]);
end

figure;errorbar(1:length(data.coincidence_indx_tk_all),...
    cellfun(@nanmean,data.coincidence_indx_tk_all),...
    ci(1,:),ci(2,:));
%%
for i = 1:length(data.coincidence_indx_all)
    figure;
    hist(data.coincidence_indx_all{i},100)
end


%%
figure;
animal = 34;
listing_cell = {...
    '111101' ...
    '111102' ...
    '111103' ...
    '111105' ...
    '111106' ...
    '111107' ...
    '111108' ...
    '111109' ...
    '111110' ...
    }'

ci = zeros(2,length(data.coincidence_indx_all));
ci_gauss = zeros(1,length(data.coincidence_indx_all));

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal data spike_num ci ci_gauss
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = [];
    load([curr_loop_folder filesep 'Manual_ROI.mat'],'roi_trace_long');
    numframes = size(roi_trace_long,2);
    for i = 1:length(data.spike_frames)
        spike_num(trace_loop_file)=sum(cellfun(@length,data.spike_frames{trace_loop_file}))/(numframes*size(roi_trace_long,1));
    end
    
    num_rep = 1000;
    ci(:,trace_loop_file) = bootci(num_rep,@(x) sum(cellfun(@length,x))/(numframes*size(roi_trace_long,1)),data.spike_frames{trace_loop_file});
    
    % gaussian error bars
    ci_gauss(trace_loop_file) = 1.96*(std(cellfun(@length,data.spike_frames{trace_loop_file}))/(numframes*size(roi_trace_long,1))...
        /sqrt(length(data.spike_frames{trace_loop_file})));
end
plot(spike_num);
errorbar(1:length(spike_num),spike_num,ci(1,:),ci(2,:));
hold on
errorbar(spike_num,ci_gauss,'r')


%% Plot the coincidence index over time
mean_coincidence = cellfun(@nanmean,data.coincidence_indx_all);
std_coincidence = cellfun(@nanstd,data.coincidence_indx_all);
ci_gauss = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x))),data.coincidence_indx_all);

figure;
errorbar(mean_coincidence,ci_gauss);
xlabel('Session')
ylabel('AP coincidence index')

mean_coincidence_tk = cellfun(@nanmean,data.coincidence_indx_tk_all);
ci_tk_gauss = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x))),data.coincidence_indx_tk_all);

figure;
errorbar(mean_coincidence_tk,ci_tk_gauss);
xlabel('Session')
ylabel('TK coincidence index')

% plot the coincidence per cell per spike - DOESN'T ACTUALLY MAKE SENSE

% num_spikes = [];
% for i = 1:length(data.spike_frames)
%     num_spikes(i) = sum(cellfun(@length,data.spike_frames{i}));
% end
% 
% coincidence_indx_all_spikenorm = {};
% for i = 1:length(data.coincidence_indx_all)
%     coincidence_indx_all_spikenorm{i} = data.coincidence_indx_all{i}/num_spikes(i);
% end
% mean_coincidence_spikenorm = cellfun(@nanmean,coincidence_indx_all_spikenorm);
% ci_spikenorm_gauss = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x))),coincidence_indx_all_spikenorm);
% figure;
% errorbar(mean_coincidence_spikenorm,ci_spikenorm_gauss);
% xlabel('Session')
% ylabel('AP coincidence index/# total spikes')
% 
% coincidence_indx_tk_all_spikenorm = {};
% for i = 1:length(data.coincidence_indx_all)
%     coincidence_indx_tk_all_spikenorm{i} = data.coincidence_indx_tk_all{i}/num_spikes(i);
% end
% mean_coincidence_tk_spikenorm = cellfun(@nanmean,coincidence_indx_tk_all_spikenorm);
% ci_tk_spikenorm_gauss = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x))),coincidence_indx_tk_all_spikenorm);
% figure;
% errorbar(mean_coincidence_tk_spikenorm,ci_tk_spikenorm_gauss);
% xlabel('Session')
% ylabel('TK coincidence index/# total spikes')

%% plot cells but with TK index

for curr_session = 1:length(data.coincidence_indx_tk_all);
    % plot the ROIs
    figure
    hold on;
    set(gca,'color','k');
    set(gcf,'color','k');
    for i = 1:length(data.polygon{curr_session}.ROI)
        r = patch(data.polygon{curr_session}.ROI{i}(:,1),-data.polygon{curr_session}.ROI{i}(:,2),[0.5 0.5 0.5]);
        %set(r,'linewidth',length(data.spike_frames{curr_session}{i})/10+0.01);
        set(r,'edgecolor','w');
    end
    
%     for i = 1:length(data.mod_cells_early{curr_session})
%         curr_color = data.mod_cells_early_percent{curr_session}(i)*0.5+0.5;
%         r = patch(data.polygon{curr_session}.ROI{data.mod_cells_early{curr_session}(i)}(:,1),-data.polygon{curr_session}.ROI{data.mod_cells_early{curr_session}(i)}(:,2),[0 curr_color 0]);
%         set(r,'linewidth',length(data.spike_frames{curr_session}{data.mod_cells_early{curr_session}(i)})/10+0.01);
%         set(r,'edgecolor','w');
%     end
%     
%     for i = 1:length(data.mod_cells_middle{curr_session})
%         curr_color = data.mod_cells_middle_percent{curr_session}(i)*0.5+0.5;
%         r = patch(data.polygon{curr_session}.ROI{data.mod_cells_middle{curr_session}(i)}(:,1),-data.polygon{curr_session}.ROI{data.mod_cells_middle{curr_session}(i)}(:,2),[curr_color curr_color 0]);
%         set(r,'linewidth',length(data.spike_frames{curr_session}{data.mod_cells_middle{curr_session}(i)})/10+0.01);
%         set(r,'edgecolor','w');
%     end
%     
%     for i = 1:length(data.mod_cells_late{curr_session})
%         curr_color = data.mod_cells_late_percent{curr_session}(i)*0.5+0.5;
%         r = patch(data.polygon{curr_session}.ROI{data.mod_cells_late{curr_session}(i)}(:,1),-data.polygon{curr_session}.ROI{data.mod_cells_late{curr_session}(i)}(:,2),[curr_color  0]);
%         set(r,'linewidth',length(data.spike_frames{curr_session}{data.mod_cells_late{curr_session}(i)})/10+0.01);
%         set(r,'edgecolor','w');
%     end
    
    ylim([-512 0]);
    xlim([0 512]);
    
    % plot lines on ROIs - thicker = more correlated
    for i = 1:length(data.coincidence_indx_all{curr_session})
        cell1 = data.cell_pair_all{curr_session}(i,1);
        cell2 = data.cell_pair_all{curr_session}(i,2);
        if data.coincidence_indx_tk_all{curr_session}(i) > 0
            center_x_1 = mean(data.polygon{curr_session}.ROI{cell1}(:,1));
            center_y_1 = mean(data.polygon{curr_session}.ROI{cell1}(:,2));
            center_x_2 = mean(data.polygon{curr_session}.ROI{cell2}(:,1));
            center_y_2 = mean(data.polygon{curr_session}.ROI{cell2}(:,2));
            x = [center_x_1 center_x_2];
            y = [center_y_1 center_y_2];
            z = [0 0];
            r = patch(x,-y,z);
            set(r,'edgealpha',data.coincidence_indx_tk_all{curr_session}(i));
            set(r,'linewidth',data.num_coincidence_all{curr_session}(i)/10);
            set(r,'edgecolor','w');
        end
    end
end

%% plot correlation for mod cells and nonmod cells
% need to flip cell orientation
mod_cells_early_flip = cellfun(@transpose,data.mod_cells_early,'UniformOutput',0);
mod_cells_middle_flip = cellfun(@transpose,data.mod_cells_middle,'UniformOutput',0);
mod_cells_late_flip = cellfun(@transpose,data.mod_cells_late,'UniformOutput',0);
mod_cells = unique([mod_cells_early_flip{:} mod_cells_middle_flip{:} mod_cells_late_flip{:}]);
% get modulated cells
mod_indx = {};
nonmod_indx = {};
for i = 1:length(data.cell_pair_all)
    mod_indx{i} = ismember(data.cell_pair_all{i}(:,1),mod_cells) & ...
        ismember(data.cell_pair_all{i}(:,2),mod_cells);
    nonmod_indx{i} = ~mod_indx{i};
end
% build coincidence indicies for mod cells
coincidence_indx_all_mod = {};
for i = 1:length(data.cell_pair_all)
    coincidence_indx_all_mod{i} = data.coincidence_indx_all{i}(mod_indx{i});
%     coincidence_indx_all_mod{i}(isnan(coincidence_indx_all_mod{i})) = 0;
end
coincidence_indx_tk_all_mod = {};
for i = 1:length(data.cell_pair_all)
    coincidence_indx_tk_all_mod{i} = data.coincidence_indx_tk_all{i}(mod_indx{i});
%     coincidence_indx_tk_all_mod{i}(isnan(coincidence_indx_tk_all_mod{i})) = 0;
end
% build coincidence indicies for nonmod cells
coincidence_indx_all_nonmod = {};
for i = 1:length(data.cell_pair_all)
    coincidence_indx_all_nonmod{i} = data.coincidence_indx_all{i}(nonmod_indx{i});
%     coincidence_indx_all_nonmod{i}(isnan(coincidence_indx_all_nonmod{i})) = 0;
end
coincidence_indx_tk_all_nonmod = {};
for i = 1:length(data.cell_pair_all)
    coincidence_indx_tk_all_nonmod{i} = data.coincidence_indx_tk_all{i}(nonmod_indx{i});
%     coincidence_indx_tk_all_nonmod{i}(isnan(coincidence_indx_tk_all_nonmod{i})) = 0;
end

% plot raw AP values
figure; hold on;
for i = 1:length(coincidence_indx_all_nonmod)
    plot(i,coincidence_indx_all_nonmod{i},'.b')
end
for i = 1:length(coincidence_indx_all_mod)
    plot(i,coincidence_indx_all_mod{i},'.r')
end
xlabel('Session')
ylabel('AP coincidence index')
legend({'NonMod','Mod'})
% plot raw TK values
figure; hold on;
for i = 1:length(coincidence_indx_tk_all_nonmod)
    plot(i,coincidence_indx_tk_all_nonmod{i},'.b')
end
for i = 1:length(coincidence_indx_tk_all_mod)
    plot(i+0.5,coincidence_indx_tk_all_mod{i},'.r')
end
xlabel('Session')
ylabel('TK coincidence index')
legend({'NonMod','Mod'})

% get mean/ci mod
mean_coincidence_mod = cellfun(@nanmean,coincidence_indx_all_mod);
ci_gauss_mod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),coincidence_indx_all_mod);
mean_coincidence_tk_mod = cellfun(@nanmean,coincidence_indx_tk_all_mod);
ci_gauss_tk_mod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),coincidence_indx_tk_all_mod);
% get mean/ci nonmod
mean_coincidence_nonmod = cellfun(@nanmean,coincidence_indx_all_nonmod);
ci_gauss_nonmod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),coincidence_indx_all_nonmod);
mean_coincidence_tk_nonmod = cellfun(@nanmean,coincidence_indx_tk_all_nonmod);
ci_gauss_tk_nonmod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),coincidence_indx_tk_all_nonmod);

% plot for AP
figure; hold on;
errorbar(mean_coincidence_nonmod,ci_gauss_nonmod);
errorbar(mean_coincidence_mod,ci_gauss_mod,'r');
xlabel('Session')
ylabel('AP coincidence index')
legend({'NonMod','Mod'})
% plot for TK
figure; hold on;
errorbar(mean_coincidence_tk_nonmod,ci_gauss_tk_nonmod);
errorbar(mean_coincidence_tk_mod,ci_gauss_tk_mod,'r');
xlabel('Session')
ylabel('TK coincidence index')
legend({'NonMod','Mod'})

% plot it with boxplot (doesn't work yet, but this is the idea)
% boxprep = [coincidence_indx_tk_all_mod{1};coincidence_indx_tk_all_nonmod{1}];
% boxgroup = [ones(length(coincidence_indx_tk_all_mod{1}),1);2*ones(length(coincidence_indx_tk_all_nonmod{1}),1)];

%% plot event rate
animal = 35;
listing_cell = {...
    '111101' ...
    }';

spike_rate_mod = {};
spike_rate_nonmod = {};

% need to flip cell orientation
mod_cells_early_flip = cellfun(@transpose,data.mod_cells_early,'UniformOutput',0);
mod_cells_middle_flip = cellfun(@transpose,data.mod_cells_middle,'UniformOutput',0);
mod_cells_late_flip = cellfun(@transpose,data.mod_cells_late,'UniformOutput',0);
mod_cells = unique([mod_cells_early_flip{:} mod_cells_middle_flip{:} mod_cells_late_flip{:}]);
nonmod_cells = setdiff([1:length(data.spike_frames{1})],mod_cells);

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal data spike_rate_mod spike_rate_nonmod mod_cells nonmod_cells
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\Data\AP' num2str(animal) '\' listing_cell{trace_loop_file}];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = [];
    load([curr_loop_folder filesep 'Manual_ROI.mat'],'roi_trace_long');
    numframes = size(roi_trace_long,2);
    for i = 1:length(data.spike_frames)
        spike_rate_mod{trace_loop_file}=sum(cellfun(@length,data.spike_frames{trace_loop_file}(mod_cells)))/(numframes*length(mod_cells));
        spike_rate_nonmod{trace_loop_file}=sum(cellfun(@length,data.spike_frames{trace_loop_file}(nonmod_cells)))/(numframes*length(nonmod_cells));
    end
end

% get mean/ci mod/nonmod
mean_spike_rate_mod = cellfun(@nanmean,spike_rate_mod);
ci_gauss_spike_rate_mod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),spike_rate_mod);
mean_spike_rate_nonmod = cellfun(@nanmean,spike_rate_nonmod);
ci_gauss_spike_rate_nonmod = cellfun(@(x) 1.96*(nanstd(x)/sqrt(length(x)/2)),spike_rate_nonmod);

figure; hold on;
errorbar(mean_spike_rate_nonmod,ci_gauss_spike_rate_nonmod);
errorbar(mean_spike_rate_mod,ci_gauss_spike_rate_mod,'r');
xlabel('Session')
ylabel('Spike rate')
legend({'NonMod','Mod'})

