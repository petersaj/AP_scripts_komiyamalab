%% Process raw data: motion correct, make summed movies

day = '151213';
animal = 'AP161';

data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
corrected_path = ['/usr/local/lab/People/Andy/Data/corticospinal_fastz_imaging/' day '_' animal];

data_path_dir = dir(data_path);
regions = {data_path_dir(3:end).name};

for curr_region = 1:length(regions);
    
    %%% Motion correct region z-stacks
    
    zstack_path = [data_path filesep regions{curr_region} filesep 'zstack'];
    
    AP_mcstack(zstack_path);
    
    % Move corrected zstack to People folder
    mc_path = [zstack_path filesep 'motion_corrected_avg'];
    corrected_path_region = [corrected_path filesep regions{curr_region}];
    mkdir(corrected_path_region);
    movefile([mc_path filesep '*'],corrected_path_region);
    
    %%% Motion correct fast z
    
    region_dir = dir([data_path filesep regions{curr_region}]);
    region_paths = {region_dir(3:end).name};
    non_zstack = cellfun(@(x) ~strcmp('zstack',x),region_paths);
    
    cell_paths = cellfun(@(x) [data_path filesep regions{curr_region} ...
        filesep x],region_paths(non_zstack),'uni',false);
    
    for curr_cell = 1:length(cell_paths)
        
        AP_mcfastz(cell_paths{curr_cell}, [2,8]);
        % Downsample motion corrected fast-z
        mc_path = [cell_paths{curr_cell} filesep 'motion_corrected'];
        AP_downsampleMovie(10,mc_path,[],true,2);
        
        % Move corrected and downsampled files to People folder
        corrected_path_cell = [corrected_path filesep regions{curr_region} filesep num2str(curr_cell)];
        mkdir(corrected_path_cell);
        movefile(mc_path,corrected_path_cell);
    end   
    
end


%% Load traces from drawn and labeled ROIs, plot and package

days = {'151212','151212'};
animals = {'AP162','AP163'};

df_all = cell(length(days),1);
cell_id_all = cell(length(days),1);

for curr_data = 1:length(days);
    
    day = days{curr_data};
    animal = animals{curr_data};
    
    corrected_path = ['/usr/local/lab/People/Andy/Data/corticospinal_fastz_imaging/' day '_' animal];
    
    data_path_dir = dir(corrected_path);
    regions = {data_path_dir(3:end).name};
    
    curr_somadendrite = 1;
    df = {};
    cell_id = {};
    
    for curr_region = 1:length(regions);
        
        region_dir = dir([corrected_path filesep regions{curr_region}]);
        region_folders = find([region_dir.isdir]);
        region_paths = {region_dir(region_folders(3:end)).name};
        cell_paths = cellfun(@(x) [corrected_path filesep regions{curr_region} ...
            filesep x],region_paths,'uni',false);
        
        for curr_cell = 1:length(cell_paths)
            
            roi_path = [cell_paths{curr_cell} filesep 'motion_corrected/summed'];
            roi_dir = dir([roi_path filesep '*.roi']);
            roi_filenames = cellfun(@(x) [roi_path filesep x],sort({roi_dir(:).name}),'uni',false);
            
            mc_path = [cell_paths{curr_cell} filesep 'motion_corrected'];
            roi_trace = AP_getConcatTrace_continuous_batch_fastz( ...
                roi_filenames,[1,2],mc_path);
            
            % These traces are short and the decay can be long, so just do
            % dirty baseline estimation
            dendrite_trace_df = cell2mat(arrayfun(@(x) (roi_trace{1}(x,:) - ...
                prctile(roi_trace{1}(x,:),10))./ ...
                prctile(roi_trace{1}(x,:),10),[1:size(roi_trace{1},1)]','uni',false));
            
            soma_trace_df = cell2mat(arrayfun(@(x) (roi_trace{2}(x,:) - ...
                prctile(roi_trace{2}(x,:),10))./ ...
                prctile(roi_trace{2}(x,:),10),[1:size(roi_trace{2},1)]','uni',false));
            
            % If multiple somas, load labels
            if size(soma_trace_df,1) > 1
                labels = load([roi_path filesep 'slice_1.roilabel'],'-MAT');
                dendrite_somas = str2num(cell2mat([labels.roi_labels{:}]'));
                % If not, all dendrites are for the only soma
            else
                dendrite_somas = ones(size(dendrite_trace_df,1),1);
            end
            
            % Store each soma/dendrite pair
            for curr_soma = 1:size(soma_trace_df,1)
                
                curr_dendrites = dendrite_somas == curr_soma;
                
                df{curr_somadendrite} = [soma_trace_df(curr_soma,:); ...
                    dendrite_trace_df(curr_dendrites,:)];
                cell_id{curr_somadendrite} = [curr_data,curr_region,curr_cell,curr_soma];
                
                curr_somadendrite = curr_somadendrite + 1;
                
            end
            
        end
    end
    
    df_all{curr_data} = df;
    cell_id_all{curr_data} = cell_id;
    
end

% Plot individual recordings
df_cat = horzcat(df_all{:});

figure;
sq = ceil(sqrt(length(df_cat)));
for curr_soma = 1:length(df_cat);
    subplot(sq,sq,curr_soma); hold on;
    curr_trace = bsxfun(@times,df_cat{curr_soma},1./max(df_cat{curr_soma},[],2));
    for curr_plot = 1:size(curr_trace,1)
        if curr_plot == 1;
            col = 'r';
        else
            col = 'k';
        end
        plot(curr_trace(curr_plot,:) + curr_plot,col);
    end
    xlim([0,size(curr_trace,2)]);
    ylim([1,size(curr_trace,1)+1]);
    axis off
end

% Save data
save('/usr/local/lab/People/Andy/Corticospinal/corticospinal_analyses/fast_z/dendrite_soma_data.mat','dendrite_soma_data');

%% Get images of dendrites and somas

% Copy over summed movies and roi files for dendrite/soma
data_path = '/usr/local/lab/People/Andy/Data/corticospinal_fastz_imaging';
animal_paths = {'151212_AP162','151212_AP163'};

cell_idx = 1;

soma_fig = figure('Name','Somas');
dendrite_fig = figure('Name','Dendrites');

for curr_animal = 1:length(animal_paths)
    
    curr_path = [data_path filesep animal_paths{curr_animal}];
    curr_dir = dir(curr_path);
    curr_regions = {curr_dir(3:end).name};
    
    for curr_region = 1:length(curr_regions)
        
        curr_region_path = [curr_path filesep curr_regions{curr_region}];
        region_dir = dir(curr_region_path);
        region_dir_folders = {region_dir([region_dir.isdir]).name};
        region_subregions = region_dir_folders(3:end);
        
        for curr_subregion = 1:length(region_subregions);
            
            curr_subregion_path = [curr_region_path filesep region_subregions{curr_subregion} filesep ...
                'motion_corrected' filesep 'summed'];
            curr_summed_filename = dir([curr_subregion_path filesep '*.tif']);
            im = AP_load_tiff([curr_subregion_path filesep curr_summed_filename.name]);
            
            curr_dendrite_im = nanmean(im(:,:,1:2:end),3);
            curr_soma_im = nanmean(im(:,:,2:2:end),3);
            
            roi_dir = dir([curr_subregion_path filesep '*.roi']);
            roi_filenames = sort({roi_dir.name});
            
            curr_dendrite_rois = load([curr_subregion_path filesep roi_filenames{1}],'-MAT');
            curr_soma_rois = load([curr_subregion_path filesep roi_filenames{2}],'-MAT');
            
            if length(curr_soma_rois.polygon.ROI) > 1;
                roilabel_filenames = dir([curr_subregion_path filesep '*.roilabel']);
                curr_roilabels = load([curr_subregion_path filesep roilabel_filenames.name],'-MAT');
                
            end
            
            for curr_cell = 1:length(curr_soma_rois.polygon.ROI)
                
                curr_soma_roi = curr_soma_rois.polygon.ROI{curr_cell};
                
                if length(curr_soma_rois.polygon.ROI) > 1;
                    curr_dendrite_roi_idx = ...
                        cellfun(@(x) strcmp(num2str(curr_cell),x),curr_roilabels.roi_labels);
                else
                    curr_dendrite_roi_idx = true(size(curr_dendrite_rois.polygon.ROI));
                end
                
                curr_dendrite_rois_cell = curr_dendrite_rois.polygon.ROI(curr_dendrite_roi_idx);
                
                % Plot plane image and ROI
                figure(dendrite_fig);
                subplot(7,7,cell_idx); hold on;
                imagesc(curr_dendrite_im);colormap(gray);
                for curr_dendrite = 1:length(curr_dendrite_rois_cell);
                    plot(curr_dendrite_rois_cell{curr_dendrite}(:,1), ...
                        curr_dendrite_rois_cell{curr_dendrite}(:,2),'b');
                end
                axis square; axis off;
                px_leeway = 20;
                curr_dendrite_cat = vertcat(curr_dendrite_rois_cell{:});
                axis_dim = ceil(max(max(curr_dendrite_cat) - ...
                    min(curr_dendrite_cat)));
                min_x = floor(min(curr_dendrite_cat(:,1)));
                min_y = floor(min(curr_dendrite_cat(:,2)));
                xlim([min_x-px_leeway,min_x+axis_dim+px_leeway]);
                ylim([min_y-px_leeway,min_y+axis_dim+px_leeway]);
                
                figure(soma_fig);
                subplot(7,7,cell_idx); hold on;
                imagesc(curr_soma_im);colormap(gray);
                plot(curr_soma_roi(:,1), ...
                    curr_soma_roi(:,2),'b');
                axis square; axis off;
                px_leeway = 20;
                axis_dim = ceil(max(max(curr_soma_roi) - ...
                    min(curr_soma_roi)));
                min_x = floor(min(curr_soma_roi(:,1)));
                min_y = floor(min(curr_soma_roi(:,2)));
                xlim([min_x-px_leeway,min_x+axis_dim+px_leeway]);
                ylim([min_y-px_leeway,min_y+axis_dim+px_leeway]);
                
                cell_idx = cell_idx + 1;
            end           
        end
    end
end 


%% Classify overlapping and unique events (by eye, unblinded)

% load data
load('/usr/local/lab/People/Andy/Corticospinal/corticospinal_analyses/fast_z/dendrite_soma_data.mat');

df_cat = horzcat(dendrite_soma_data.df{:});

dendrite_soma_events_all = cell(length(df_cat),1);
soma_unique_events_all = cell(length(df_cat),1);
alldendrite_unique_events_all = cell(length(df_cat),1);
dendrite_unique_events_all = cell(length(df_cat),1);

f = figure;
for curr_soma = 1:length(df_cat);
    
    set(f,'Name',[num2str(curr_soma) '/' num2str(length(df_cat))]);
    
    % Plot data
    clf(f);
    curr_trace = bsxfun(@times,df_cat{curr_soma},1./max(df_cat{curr_soma},[],2));
    p = nan(size(curr_trace,1),1);
    for curr_plot = 1:size(curr_trace,1)
        if curr_plot == 1
            col = 'b';
        else
            col = 'k';
        end
        p(curr_plot) = subplot(size(curr_trace,1),1,curr_plot);
        plot(curr_trace(curr_plot,:),col);
        xlim([0,size(curr_trace,2)]);
        ylim([min(curr_trace(curr_plot,:)),max(curr_trace(curr_plot,:))]);
        axis off
    end
    linkaxes(p,'x');
    
    % Manually set event times
   
    clear dendrite_soma_events soma_unique_events alldendrite_unique_events dendrite_unique_events;   
    
    % Shared dendrite/soma events
    title(p(1),'Identify shared dendrite/soma events');
    [dendrite_soma_events,~] = ginput;
    for i = 1:size(dendrite_soma_events,1)
        for curr_p = 1:length(p)
            subplot(p(curr_p));
            line(repmat(dendrite_soma_events(i),1,2),ylim,'color','r')
        end        
    end
    
    % Unique soma events
    title(p(1),'Identify soma-only events');
    [soma_unique_events,~] = ginput;
    for i = 1:size(soma_unique_events,1)
        for curr_p = 1:length(p)
            subplot(p(curr_p));
            line(repmat(soma_unique_events(i),1,2),ylim,'color','b')
        end        
    end
    
    % Unique shared dendrite events
    title(p(1),'Identify shared dendrite-only events');
    [alldendrite_unique_events,~] = ginput;
    for i = 1:size(alldendrite_unique_events,1)
        for curr_p = 1:length(p)
            subplot(p(curr_p));
            line(repmat(alldendrite_unique_events(i),1,2),ylim,'color','k')
        end        
    end
    
    % Unique dendrite events
    title(p(1),'Identify unique dendrite events');
    [dendrite_unique_events,~] = ginput;
    for i = 1:size(dendrite_unique_events,1)
        for curr_p = 1:length(p)
            subplot(p(curr_p));
            line(repmat(dendrite_unique_events(i),1,2),ylim,'color','k','linestyle','--')
        end        
    end
    
    dendrite_soma_events_all{curr_soma} = dendrite_soma_events;
    soma_unique_events_all{curr_soma} = soma_unique_events;
    alldendrite_unique_events_all{curr_soma} = alldendrite_unique_events;
    dendrite_unique_events_all{curr_soma} = dendrite_unique_events;
    
    title(p(1),'Final events');
    waitforbuttonpress;
    
end

dendrite_soma_events_manual.dendrite_soma_events = dendrite_soma_events_all;
dendrite_soma_events_manual.soma_unique_events = soma_unique_events_all;
dendrite_soma_events_manual.alldendrite_unique_events = alldendrite_unique_events_all;
dendrite_soma_events_manual.dendrite_unique_events = dendrite_unique_events_all;

% Save data
save('/usr/local/lab/People/Andy/Corticospinal/corticospinal_analyses/fast_z/dendrite_soma_events_manual.mat','dendrite_soma_events_manual');


%% Load, quantify, plot manual events

% Load data
load('/usr/local/lab/People/Andy/Corticospinal/corticospinal_analyses/fast_z/dendrite_soma_events_manual.mat');
% Load events
load('/usr/local/lab/People/Andy/Corticospinal/corticospinal_analyses/fast_z/dendrite_soma_data.mat');

% Fraction of events which are shared, soma-specific, dendrite-specific
all_events = [dendrite_soma_events_manual.dendrite_soma_events, ...
    dendrite_soma_events_manual.soma_unique_events, ...
    dendrite_soma_events_manual.alldendrite_unique_events, ...
    dendrite_soma_events_manual.dendrite_unique_events];

all_events_n = cellfun(@length,all_events);

shared_frac = all_events_n(:,1)./sum(all_events_n,2);
soma_unique_frac = all_events_n(:,2)./sum(all_events_n,2);
alldendrite_unique_frac = all_events_n(:,3)./sum(all_events_n,2);
dendrite_unique_frac = all_events_n(:,4)./sum(all_events_n,2);

figure; hold on;
cat_data = [shared_frac,soma_unique_frac,alldendrite_unique_frac,dendrite_unique_frac];
frac_mean = nanmean(cat_data,1);
frac_sem = nanstd(cat_data,[],1)./sqrt(sum(~isnan(cat_data),1));
bar(frac_mean,'FaceColor','w','linewidth',2);
errorbar(frac_mean,frac_sem,'.k','linewidth',2);

plotSpread(cat_data,0.1);

ylabel('Fraction of events');
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'Shared','Soma-specific','Dendrite-specific','Dendrite branch specific'});

% Plot all data with marked events
df_cat = horzcat(dendrite_soma_data.df{:});

figure;
sq = ceil(sqrt(length(df_cat)));
for curr_soma = 1:length(df_cat);
    
    % Plot traces
    subplot(sq,sq,curr_soma); hold on;
    curr_trace = bsxfun(@times,df_cat{curr_soma},1./max(df_cat{curr_soma},[],2));
    for curr_plot = 1:size(curr_trace,1)
        if curr_plot == 1;
            col = 'r';
        else
            col = 'k';
        end
        plot(curr_trace(curr_plot,:) + curr_plot,col);
    end
    xlim([0,size(curr_trace,2)]);
    ylim([1,size(curr_trace,1)+1]);
    axis off
    
    % Plot events
    curr_shared_events = dendrite_soma_events_manual.dendrite_soma_events{curr_soma};
    curr_soma_events = dendrite_soma_events_manual.soma_unique_events{curr_soma};
    curr_dendrite_events = dendrite_soma_events_manual.alldendrite_unique_events{curr_soma};
    curr_dendrite_branch_events = dendrite_soma_events_manual.dendrite_unique_events{curr_soma};
    
    for i = 1:length(curr_shared_events)
        line(repmat(curr_shared_events(i),1,2),ylim,'color','b');
    end
    for i = 1:length(curr_soma_events)
        line(repmat(curr_soma_events(i),1,2),ylim,'color','r');
    end
    for i = 1:length(curr_dendrite_events)
        line(repmat(curr_dendrite_events(i),1,2),ylim,'color','k');
    end
    for i = 1:length(curr_dendrite_branch_events)
        line(repmat(curr_dendrite_branch_events(i),1,2),ylim,'color','k','linestyle','--');
    end
end

% Get raw numbers for text
disp(['All events: ' num2str(sum(all_events_n,1))]);


%% Get within/across cell correlations from long movies of z-stack animals

pairwise_corr_all = cell(2,2);

% ANIMAL 1

% Set paths
data_path = '/usr/local/lab/People/Andy/Data/AP162';
roi_labels_filename = [data_path filesep 'AP162_roilabels.roilabel'];
data_filename = [data_path filesep 'AP162_rois' filesep ...
    '151209_AP162_001_001_corrected_summed_50_analysis'];

% Load data
data = load(data_filename);

% Mark bad ROIs manually
bad_rois = [17,40,12];

% Load branch labels
curr_labels = load(roi_labels_filename,'-MAT');
branched_rois = cellfun(@(x) ~isempty(x),curr_labels.roi_labels);
num_labels = nan(length(curr_labels.roi_labels),1);
num_labels(branched_rois) = ...
    cellfun(@(x) str2num(x{1}),curr_labels.roi_labels(branched_rois));
num_labels(bad_rois) = [];

pairwise_branch = AP_itril(repmat(num_labels,1,length(num_labels)) == ...
    repmat(num_labels',length(num_labels),1),-1);

% Smooth data
curr_smoothpeak = nan(size(data.im.roi_trace_df));
for curr_roi = 1:size(curr_smoothpeak,1)
    curr_smoothpeak(curr_roi,:) = smooth(data.im. ...
        roi_trace_df(curr_roi,:),100,'loess');
end

% Get normalized dot product
curr_smoothpeak_norm = ...
    bsxfun(@times,curr_smoothpeak,1./sqrt(sum(curr_smoothpeak.^2,2)));

curr_smoothpeak_norm(bad_rois,:) = [];
pairwise_corr = AP_itril(curr_smoothpeak_norm*curr_smoothpeak_norm',-1);

pairwise_corr_all(1,:) = {pairwise_corr(pairwise_branch),pairwise_corr(~pairwise_branch)};

% ANIMAL 2

% Set paths
data_path = '/usr/local/lab/People/Andy/Data/AP163';
roi_labels_filename = [data_path filesep 'AP163_roilabels.roilabel'];
data_filename = [data_path filesep 'AP163_rois' filesep ...
    '151209_AP163_001_001_corrected_summed_50_analysis'];

% Load data
data = load(data_filename);

% Mark bad ROIs manually
bad_rois = [12,31,20,43,17,27,37,38,55,47,48,35];

% Load branch labels
curr_labels = load(roi_labels_filename,'-MAT');
branched_rois = cellfun(@(x) ~isempty(x),curr_labels.roi_labels);
num_labels = nan(length(curr_labels.roi_labels),1);
num_labels(branched_rois) = ...
    cellfun(@(x) str2num(x{1}),curr_labels.roi_labels(branched_rois));
num_labels(bad_rois) = [];

pairwise_branch = AP_itril(repmat(num_labels,1,length(num_labels)) == ...
    repmat(num_labels',length(num_labels),1),-1);

% Smooth data
curr_smoothpeak = nan(size(data.im.roi_trace_df));
for curr_roi = 1:size(curr_smoothpeak,1)
    curr_smoothpeak(curr_roi,:) = smooth(data.im. ...
        roi_trace_df(curr_roi,:),100,'loess');
end

% Get normalized dot product
curr_smoothpeak_norm = ...
    bsxfun(@times,curr_smoothpeak,1./sqrt(sum(curr_smoothpeak.^2,2)));

curr_smoothpeak_norm(bad_rois,:) = [];
pairwise_corr = AP_itril(curr_smoothpeak_norm*curr_smoothpeak_norm',-1);

pairwise_corr_all(2,:) = {pairwise_corr(pairwise_branch),pairwise_corr(~pairwise_branch)};

% Plot correlation of all pairs
figure;

subplot(2,2,1);
plotSpread(pairwise_corr_all(1,:));
line(xlim,[0.8,0.8],'color','r','linestyle','--');
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Sibling','Non-sibling'});
ylim([-1 1]);
ylabel('Normalized dot product');
title('AP162');

subplot(2,2,2);
plotSpread(pairwise_corr_all(2,:));
line(xlim,[0.8,0.8],'color','r','linestyle','--');
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Sibling','Non-sibling'});
ylim([-1 1]);
ylabel('Normalized dot product');
title('AP163');

grp_sib = vertcat(pairwise_corr_all{:,1});
grp_nonsib = vertcat(pairwise_corr_all{:,2});

subplot(2,2,3);
plotSpread({grp_sib,grp_nonsib});
line(xlim,[0.8,0.8],'color','r','linestyle','--');
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Sibling','Non-sibling'});
ylim([-1 1]);
ylabel('Normalized dot product');
title('Combined');

subplot(2,2,4);
boxplot([grp_sib;grp_nonsib],[ones(size(grp_sib));2*ones(size(grp_nonsib))]);
line(xlim,[0.8,0.8],'color','r','linestyle','--');
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Sibling','Non-sibling'});
ylim([-1 1]);
ylabel('Normalized dot product');
title('Combined');

% Plot histograms for each population indivually normalized
bin_edges = linspace(-1,1,100);
bin_centers = [bin_edges(1:end-1) + diff(bin_edges)/2];
bin_edges(end) = Inf;

grp_sib_hist = histc(grp_sib,bin_edges)./length(grp_sib);
grp_nonsib_hist = histc(grp_nonsib,bin_edges)./length(grp_nonsib);

figure; 

subplot(2,1,1);hold on;
bar(bin_centers,grp_nonsib_hist(1:end-1),'FaceColor','k','EdgeColor','none','BarWidth',1);
bar(bin_centers,grp_sib_hist(1:end-1),'FaceColor','r','EdgeColor','none','BarWidth',1);
xlabel('Normalized dot product')
ylabel('Fraction of pairs')
line([0.8,0.8],ylim,'color','k','linestyle','--','linewidth',2);

subplot(2,1,2);hold on;
bar(bin_centers,grp_nonsib_hist(1:end-1)./max(grp_nonsib_hist), ...
    'FaceColor','k','EdgeColor','none','BarWidth',1);
bar(bin_centers,grp_sib_hist(1:end-1)./max(grp_sib_hist), ...
    'FaceColor','r','EdgeColor','none','BarWidth',1);
xlabel('Normalized dot product')
ylabel('Normalized fraction of pairs')
line([0.8,0.8],ylim,'color','k','linestyle','--','linewidth',2);



%% For comparison: load in normal animal data, get grouping dot prod

clear all

animals = {'AP140' 'AP141' 'AP142' 'AP145' 'AP146' 'AP147' 'AP148' 'AP150'};

all_roi_dot = cell(length(animals),1);

% Grab anaylsis files (assume dendrites in specific folder)
all_analysis_path = cell(size(animals));
all_analysis_files = cell(size(animals));
if ~ispc
    data_dir = '/usr/local/lab/People/Andy/Data';
else
    data_dir = 'C:\Users\Andy\Documents\KomiyamaLab\data\analysis_files';
end
for curr_animal = 1:length(animals)
    if ~ispc
        curr_analysis_dir = [data_dir filesep animals{curr_animal} filesep ...
            animals{curr_animal} '_batch_thresh_roi'];
    else
        curr_analysis_dir = [data_dir filesep animals{curr_animal}];
    end
    
    curr_analysis_files = dir([curr_analysis_dir filesep '*.mat']);
    
    all_analysis_path{curr_animal} = curr_analysis_dir;
    all_analysis_files{curr_animal} = sort({curr_analysis_files.name});   
end

for curr_animal = 1:length(animals)
    
    animal = animals{curr_animal};
    disp(['Loading animal ' animal]);
    
    data(curr_animal).animal = animal;
    
    analysis_path = all_analysis_path{curr_animal};
    analysis_files = sort(all_analysis_files{curr_animal});
    num_sessions = length(analysis_files);
    
    
    % Load ROI data
    
    disp('Loading data...')
    
    for curr_session = 1:num_sessions
        
        curr_file = [analysis_path filesep analysis_files{curr_session}];
        curr_data = load(curr_file);
        
        % Store data
        data(curr_animal).im(curr_session) = curr_data.im;
        
    end
    
    % Get bad ROIs from manual quality control, discard
    
    field_qc_path = '/usr/local/lab/People/Andy/Data/field_qc';
    field_qc_dir = dir([field_qc_path filesep '*.mat']);
    
    curr_animal_field_qc = find(cellfun(@(x) strncmp(animal,x,5),{field_qc_dir.name}));
    if ~isempty(curr_animal_field_qc)
        
        curr_bad_rois = load([field_qc_path filesep ...
            field_qc_dir(curr_animal_field_qc).name]);
        
        for curr_session = 1:length(data(curr_animal).im)
            data(curr_animal).im(curr_session).roi_trace(curr_bad_rois.bad_rois,:) = [];
            data(curr_animal).im(curr_session).roi_trace_bg(curr_bad_rois.bad_rois,:) = [];
            data(curr_animal).im(curr_session).roi_trace_df(curr_bad_rois.bad_rois,:) = [];
        end
    end
    data(curr_animal).roi_idx = find(~curr_bad_rois.bad_rois);
    
    
    % Find likely same dendrite ROIs / consistent high correlation
        
    disp('Getting ROI correlations...')
    
    n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
    
    % Loop through days, get L2 normalization value (avoid concat data)
    cat_sum = nan(n_rois,num_sessions);
    for curr_day = 1:num_sessions
        
        curr_smoothpeak = nan(size(data(curr_animal).im(curr_day).roi_trace_df));
        for curr_roi = 1:size(curr_smoothpeak,1)
            curr_smoothpeak(curr_roi,:) = smooth(data(curr_animal).im(curr_day). ...
                roi_trace_df(curr_roi,:),100,'loess');
        end
        cat_sum(:,curr_day) = sum(curr_smoothpeak.^2,2);
    end
    % Loop through days, get total correlation between ROIs
    l2norm = sqrt(sum(cat_sum,2));
    roi_corr = zeros(n_rois,n_rois);
    for curr_day = 1:num_sessions
        
        curr_smoothpeak = nan(size(data(curr_animal).im(curr_day).roi_trace_df));
        for curr_roi = 1:size(curr_smoothpeak,1)
            curr_smoothpeak(curr_roi,:) = smooth(data(curr_animal).im(curr_day). ...
                roi_trace_df(curr_roi,:),100,'loess');
        end
        
        curr_smoothpeak = bsxfun(@times,curr_smoothpeak,1./l2norm);
        roi_corr = roi_corr + curr_smoothpeak*curr_smoothpeak';
    end
        
    all_roi_dot{curr_animal} = AP_itril(roi_corr,-1);
    
end

hi_thresh = 0.4;
all_roi_dot_overthresh = cellfun(@(x) x(x > hi_thresh),all_roi_dot,'uni',false);

all_roi_dot_grp = cellfun(@(x,y) ones(size(x))*y, ...
    all_roi_dot,num2cell(1:length(all_roi_dot))','uni',false); 

figure;
plotSpread(all_roi_dot_overthresh);
line(xlim,[0.8,0.8],'color','r','linestyle','--');
ylim([hi_thresh,1]);
xlabel('Animal');
ylabel('Normalized dot product');
title('ROI grouping cutoff');



