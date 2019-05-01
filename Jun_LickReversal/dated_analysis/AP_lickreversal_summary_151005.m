%% Correlation of average activity in cont1/2 early/late

% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

all_active_trials = cell(length(animals),1);
all_contingencies = cell(length(animals),1);
all_latereversal = cell(length(animals),1);

for curr_animal = 1:length(mice)
        
    for curr_session = 1:length(data_all(curr_animal).im);
        
        % define reversal (use only the first so that it's settled in)
        for curr_rev = 1:length(analysis.epoch_contingencies{curr_animal}{curr_session});
            
            revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
                size(analysis.odor_trials{curr_animal}{curr_session},1)];
            % if there are < 10 trials in reversal, assume early/not real
            % reversal and use the next one
            rev_trials = diff(revs);
            if rev_trials(curr_rev) < 10
                curr_rev = curr_rev + 1;
            end
            
            % only use first 2 seconds (odor delivery) of activity
            odor_frames = data_all(curr_animal).im(curr_session).framerate*2;
            curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
                analysis.epoch_frames{curr_animal}{curr_session}(2);
            curr_odor_frames = curr_frames > 0 & curr_frames <= odor_frames;
            
            curr_act = analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
                (:,curr_odor_frames,:);
            % skip this if less than 20 trials
            if size(curr_act,1) < 20
                continue
            end
            
            curr_active_trials = permute(any(curr_act,2),[1 3 2]);
            
            curr_contingency = analysis.epoch_contingencies{curr_animal}{curr_session}(curr_rev);           
            % if control animal, don't use anything past session 10 or if
            % contingency 2 has appeared (because reversal in these
            % animals, but it's late - I'm not using this data for now)
            if ismember(curr_animal,find(ctrl_animals));
                if curr_session > 10 || curr_contingency == 2
                    continue
                end
            end
            
            all_active_trials{curr_animal} = [all_active_trials{curr_animal};curr_active_trials];
            all_contingencies{curr_animal} = [all_contingencies{curr_animal}; ...
                repmat(curr_contingency, size(curr_active_trials,1),1)];
            
            % 'Late' for reversal animals = anything including and after
            % the second reversal
            if ismember(curr_animal,find(rev_animals));
                all_latereversal{curr_animal} = [all_latereversal{curr_animal}; ...
                    repmat(sum(diff(all_contingencies{curr_animal}) ~= 0) > 1, ...
                    size(curr_active_trials,1),1)];
            else
                % 'Late for control animals = anything after the 6th session
                all_latereversal{curr_animal} = [all_latereversal{curr_animal}; ...
                    repmat(curr_session >= 6, ...
                    size(curr_active_trials,1),1)];
            end
            
            
                       
        end
    end
end

cont_earlylate_corr = nan(4,4,25);
for curr_animal = find(rev_animals);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    curr_cont = all_contingencies{curr_animal};
    curr_late = all_latereversal{curr_animal};
    
    cont_earlylate_corr(:,:,curr_animal) = corrcoef([ ...
        nanmean(curr_act(curr_cont == 1 & ~curr_late,:),1);
        nanmean(curr_act(curr_cont == 2 & ~curr_late,:),1);
        nanmean(curr_act(curr_cont == 1 & curr_late,:),1);
        nanmean(curr_act(curr_cont == 2 & curr_late,:),1);]');
    
end
for curr_animal = find(ctrl_animals);
    
    curr_act = horzcat(all_active_trials{curr_animal,:});
    
    curr_late = logical(all_latereversal{curr_animal});
    
    cont_earlylate_corr(1:2,1:2,curr_animal) = corrcoef([ ...
        nanmean(curr_act(curr_late,:),1);
        nanmean(curr_act(~curr_late,:),1)]');
    
end

figure;
subplot(2,3,1);
imagesc(nanmean(cont_earlylate_corr(:,:,rev_animals & pmm_animals),3));colormap(gray)
title('PMM rev');
set(gca,'XTickLabel',{'Early 1','Early 2','Late 1','Late 2'})
set(gca,'YTickLabel',{'Early 1','Early 2','Late 1','Late 2'})
subplot(2,3,2);
imagesc(nanmean(cont_earlylate_corr(:,:,rev_animals & alm_animals),3));colormap(gray)
title('ALM rev');
set(gca,'XTickLabel',{'Early 1','Early 2','Late 1','Late 2'})
set(gca,'YTickLabel',{'Early 1','Early 2','Late 1','Late 2'})
subplot(2,3,4);
imagesc(nanmean(cont_earlylate_corr(:,:,ctrl_animals & pmm_animals),3));colormap(gray)
title('PMM ctrl');
set(gca,'XTickLabel',{'Early 1a','Early 1b','',''})
set(gca,'YTickLabel',{'Early 1a','Early 1b','',''})
subplot(2,3,5);
imagesc(nanmean(cont_earlylate_corr(:,:,ctrl_animals & alm_animals),3));colormap(gray)
title('PMM ctrl');
set(gca,'XTickLabel',{'Early 1a','Early 1b','',''})
set(gca,'YTickLabel',{'Early 1a','Early 1b','',''})

earlychange_pmm = [nanmean(cont_earlylate_corr(1,2,rev_animals & pmm_animals),3), ...
    nanmean(cont_earlylate_corr(1,2,ctrl_animals & pmm_animals),3)];
earlychange_alm = [nanmean(cont_earlylate_corr(1,2,rev_animals & alm_animals),3)
    nanmean(cont_earlylate_corr(1,2,ctrl_animals & alm_animals),3)];

subplot(2,3,[3,6]); hold on;
plot(earlychange_pmm,'k','linewidth',2);
plot(earlychange_alm,'r','linewidth',2);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Rev','Ctrl'});
xlim([0 3]);
ylabel('Correlation');
title('Early change');
legend({'PMM','ALM'});



%% Average activity timing relative to lick onset

%% Slide 6) Activity aligned/sorted by lick onset


% Get animal group indicies
rev_animals = cellfun(@(x) strcmp(x,'y'),{mice.switch});
ctrl_animals = cellfun(@(x) strcmp(x,'n'),{mice.switch});
pmm_animals = cellfun(@(x) strcmp(x,'PMM'),{mice.location});
alm_animals = cellfun(@(x) strcmp(x,'ALM'),{mice.location});

lick_aligned_maxtime = cell(length(data_all),1);

for curr_animal = 1:length(mice)
    
    lick_aligned_maxtime{curr_animal} = ...
        nan(size(analysis.epoch_activity_CL{curr_animal}{2}{1},3), ...
        length(data_all(curr_animal).im));
        
    for curr_session = 1:length(data_all(curr_animal).im)
        
        % define reversal
        curr_rev = 1;
        revs = [1;find(diff(analysis.odor_trials{curr_animal}{curr_session}(:,2)))+1;...
            size(analysis.odor_trials{curr_animal}{curr_session},1)];
        % if there are < 10 trials in reversal, assume early/not real
        % reversal and use the next one
        rev_trials = diff(revs);
        if rev_trials(curr_rev) < 10
            curr_rev = curr_rev + 1;
        end
        
        % lick isn't segregated by reversal, do that here  
        curr_trials = false(size(analysis.odor_trials{curr_animal}{curr_session},1),1);
        curr_trials(revs(curr_rev):revs(curr_rev+1)-1) = true;
        
        curr_lick = analysis.lick_rate{curr_animal}{curr_session}( ...
            curr_trials & analysis.condition_trials{curr_animal}{curr_session}(:,1),:);
        
        % Get the first time where licking is over 3 hz, sort
        [~,lick_start] = max(curr_lick > 3,[],2);
        [~,lick_start_sort] = sort(lick_start);
        
        % Get rid of weird outlier lick trials (for now just cutoff)
        early_cutoff = 50; % depends on lick rate resolution, now 1/100 sec (*10 = ms)
        early_lick_trials = lick_start < early_cutoff;
               
        % aligning activity to lick onset (v odor onset)
        
        % how to align pre/post 
        post_frames = round(data_all(curr_animal).im(curr_session).framerate*1);
        pre_frames = round(data_all(curr_animal).im(curr_session).framerate*1);
        
        % convert lick times into frames
        lick_start_frames = round(lick_start/100* ...
            data_all(curr_animal).im(curr_session).framerate);
        
        curr_frames = analysis.epoch_frames{curr_animal}{curr_session}(1) : ...
            analysis.epoch_frames{curr_animal}{curr_session}(2);
        curr_start_frame = find(curr_frames > 0,1);
        
        act_use_frames = arrayfun(@(x) (curr_start_frame + lick_start_frames(x) - pre_frames): ...
            (curr_start_frame + lick_start_frames(x) + post_frames),1:length(lick_start_frames),'uni',false);
        
        curr_act_lickalign = cell2mat( ...
            arrayfun(@(x) analysis.epoch_activity_CL{curr_animal}{curr_session}{curr_rev} ...
            (x,act_use_frames{x},:),find(~early_lick_trials),'uni',false));
        
        curr_act_lickalign_mean = permute(nanmean(curr_act_lickalign,1),[3,2,1]);
        
        [max_act,max_act_frame] = max(curr_act_lickalign_mean,[],2);
        
        max_act_frame(max_act == 0) = NaN;
        max_act_frame = max_act_frame - pre_frames;
        max_act_time = max_act_frame./data_all(curr_animal).im(curr_session).framerate;
        
        if ~isempty(max_act_time)
            lick_aligned_maxtime{curr_animal}(:,curr_session) = max_act_time;
        end

    end
    
        disp(curr_animal);
    
end


% Get significantly active cells from baseline

framerate_animals =  arrayfun(@(x) {data_all(x).im(:).framerate},1:length(data_all),'uni',false);

epoch_seconds = cellfun(@(x,y) cellfun(@(a,b) (a(1):a(2))/b,x,y,'uni',false), ...
    analysis.epoch_frames,framerate_animals,'uni',false);

baseline_sigcells = ...
    cellfun(@(odor,seconds) cellfun(@(odor,seconds) cellfun(@(odor) ...
    AP_ranksum_matrix( ...
    permute(nanmean(odor(:,seconds < 0,:),2),[1 3 2]), ...
    permute(nanmean(odor(:,seconds >= 0 & ...
    seconds < 2,:),2),[1 3 2])) ...
    < 0.05, odor,'uni',false),odor,seconds,'uni',false), ...
    analysis.epoch_activity_CL,epoch_seconds,...
    'uni',false);

% Plot
figure;

subplot(2,2,1); hold on;
curr_animals = find(rev_animals & pmm_animals);
col = lines(max(curr_animals));
for curr_animal = curr_animals
    for curr_session = 1:size(lick_aligned_maxtime{curr_animal},2)
        curr_cells = any(vertcat(baseline_sigcells{curr_animal}{curr_session}{:}),1);
        if(any(curr_cells))
            plot(curr_session, lick_aligned_maxtime{curr_animal}(curr_cells,curr_session), ...
                '.','color',col(curr_animal,:));
        end
    end
end
xlabel('Session');
ylabel('Time from lick onset to max mean act')

subplot(2,2,2); hold on;
curr_animals = find(ctrl_animals & pmm_animals);
col = lines(max(curr_animals));
for curr_animal = curr_animals
    for curr_session = 1:size(lick_aligned_maxtime{curr_animal},2)
        curr_cells = any(vertcat(baseline_sigcells{curr_animal}{curr_session}{:}),1);
        if(any(curr_cells))
            plot(curr_session, lick_aligned_maxtime{curr_animal}(curr_cells,curr_session), ...
                '.','color',col(curr_animal,:));
        end
    end
end
xlabel('Session');
ylabel('Time from lick onset to max mean act')

subplot(2,2,3); hold on;
curr_animals = find(rev_animals & alm_animals);
col = lines(max(curr_animals));
for curr_animal = curr_animals
    for curr_session = 1:size(lick_aligned_maxtime{curr_animal},2)
        curr_cells = any(vertcat(baseline_sigcells{curr_animal}{curr_session}{:}),1);
        if(any(curr_cells))
            plot(curr_session, lick_aligned_maxtime{curr_animal}(curr_cells,curr_session), ...
                '.','color',col(curr_animal,:));
        end
    end
end
xlabel('Session');
ylabel('Time from lick onset to max mean act')

subplot(2,2,4); hold on;
curr_animals = find(ctrl_animals & alm_animals);
col = lines(max(curr_animals));
for curr_animal = curr_animals
    for curr_session = 1:size(lick_aligned_maxtime{curr_animal},2)
        curr_cells = any(vertcat(baseline_sigcells{curr_animal}{curr_session}{:}),1);
        if(any(curr_cells))
            plot(curr_session, lick_aligned_maxtime{curr_animal}(curr_cells,curr_session), ...
                '.','color',col(curr_animal,:));
        end
    end
end
xlabel('Session');
ylabel('Time from lick onset to max mean act')

















