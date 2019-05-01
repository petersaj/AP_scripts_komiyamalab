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

























