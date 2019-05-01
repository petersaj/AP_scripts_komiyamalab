%% Plot discrimination behavior (THIS IS WRONG! see next cell)

for curr_animal = 1:length(animals);
    
    rewarded_odor = cell(length(data_all(curr_animal).im),1);
    contingency_reversals = cell(length(data_all(curr_animal).im),1);
    bhv_matrix = cell(length(data_all(curr_animal).im),1);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        
        contingency_reversals{curr_session} = ...
            [0;diff(data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials))] ~= 0;
        
        bhv_matrix{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end
    
    session_num_trials = cellfun(@length,contingency_reversals);

    % Grab seperately: hit rate CL/(CL+IR), false alarm rate IL/(IL+CR)
    trial_range = 10; % surrounding trials to estimate by
    trial_conv = ones(trial_range,1);
    bhv_matrix_conv = cellfun(@(x) conv2(+x,trial_conv,'same'),bhv_matrix,'uni',false);
    hit_rate = cellfun(@(x) x(:,1)./(x(:,1)+x(:,4)),bhv_matrix_conv,'uni',false);
    false_alarm_rate = cellfun(@(x) x(:,3)./(x(:,3)+x(:,2)),bhv_matrix_conv,'uni',false);
    
    % Grab d'
    
    % get interpolated rate of conditions
    bhv_matrix_rate = cell(size(bhv_matrix));
    for curr_session = 1:length(bhv_matrix);
        for condition = 1:size(bhv_matrix{curr_session},2);
            
            curr_bin = bhv_matrix{curr_session}(:,condition);
            
            % if there's none of a current condition, rate is 0
            if ~any(curr_bin)
                bhv_matrix_rate{curr_session}(:,condition) = ...
                    zeros(length(curr_bin)-1,1);
                continue
            end
            
            curr_interp = cumsum(curr_bin);
            % set all zeros (except for the first and after last occurance) to NaN;
            curr_interp(find(curr_bin(2:find(curr_bin == 1,1,'last')) == 0) + 1) = NaN;
            curr_interp(isnan(curr_interp)) = ...
                interp1(find(~isnan(curr_interp)), ...
                curr_interp(~isnan(curr_interp)), ...
                find(isnan(curr_interp)));
            
            bhv_matrix_rate{curr_session}(:,condition) = diff(curr_interp);
        end
    end
    
    % smooth the rate (only in the forward direction for quick shifts)
    rate_smooth = 20;
    rate_smooth_kernel = ones(rate_smooth,1)/rate_smooth;
    bhv_matrix_rate_smooth_complete = cellfun(@(x) conv2(x,rate_smooth_kernel),bhv_matrix_rate,'uni',false);
    bhv_matrix_rate_smooth = cellfun(@(x) x(rate_smooth-1:end,:),bhv_matrix_rate_smooth_complete,'uni',false);
    
    % d' = zscore(hit rate (correct lick/correct lick + incorrect rejection)) -
    % zscore(false alarm rate (incorrect lick / incorrect lick + correct
    % rejection))
    bhv_raw_hit_falsealarm = cellfun(@(x) [x(:,1)./(x(:,1) + x(:,4)) ...
        x(:,3)./(x(:,3) + x(:,2))],bhv_matrix_rate_smooth,'uni',false);
    bhv_all_mean = nanmean(vertcat(bhv_raw_hit_falsealarm{:}));
    bhv_all_std = nanstd(vertcat(bhv_raw_hit_falsealarm{:}));
    bhv_dprime = cellfun(@(x) diff(fliplr(bsxfun(@times,bsxfun(@minus,x, ...
        bhv_all_mean),bhv_all_std)),[],2),bhv_raw_hit_falsealarm,'uni',false);
    % this was to zscore by session
    %bhv_dprime = cellfun(@(x) diff(fliplr(bsxfun(@times,bsxfun(@minus,x, ...
    %    nanmean(x,1)),nanstd(x,[],1))),[],2),bhv_raw_hit_falsealarm,'uni',false);
    
    hit_rate_odor1cat = vertcat(hit_rate{:});
    hit_rate_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    hit_rate_odor2cat = vertcat(hit_rate{:});
    hit_rate_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    false_alarm_rate_odor1cat = vertcat(false_alarm_rate{:});
    false_alarm_rate_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    false_alarm_rate_odor2cat = vertcat(false_alarm_rate{:});
    false_alarm_rate_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    bhv_dprime_odor1cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    bhv_dprime_odor2cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    figure;
    p1 = subplot(3,1,1); hold on;
    plot(hit_rate_odor1cat,'k');
    plot(hit_rate_odor2cat,'r');
    ylabel('Hit rate')
    title(['Behavior: ' animals{curr_animal}])
    xlim([0 sum(session_num_trials)]);
    ylim([-0.2 1.2])
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color','b','linestyle','--');
    end
    p2 = subplot(3,1,2); hold on;
    plot(false_alarm_rate_odor1cat,'k');
    plot(false_alarm_rate_odor2cat,'r');
    ylabel('False alarm rate')
    xlim([0 sum(session_num_trials)]);
    ylim([-0.2 1.2])
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color','b','linestyle','--');
    end
    p3 = subplot(3,1,3); hold on;
    plot(bhv_dprime_odor1cat,'k');
    plot(bhv_dprime_odor2cat,'r');
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color','b','linestyle','--');
    end
    xlim([0 sum(session_num_trials)]);
    ylabel('D''');
    xlabel('Trial')
    
    linkaxes([p1 p2 p3],'x');
end

%% Plot discrimination behavior (with chance conf intervals line)

for curr_animal = 1:length(animals);
    
    rewarded_odor = cell(length(data_all(curr_animal).im),1);
    contingency_reversals = cell(length(data_all(curr_animal).im),1);
    bhv_matrix = cell(length(data_all(curr_animal).im),1);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        
        contingency_reversals{curr_session} = ...
            [0;diff(data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials))] ~= 0;
        
        bhv_matrix{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end
    
    session_num_trials = cellfun(@length,contingency_reversals);

    % Grab seperately: hit rate CL/(CL+IR), false alarm rate IL/(IL+CR)
    trial_range = 10; % surrounding trials to estimate by
    trial_conv = ones(trial_range,1);
    bhv_matrix_conv = cellfun(@(x) conv2(+x,trial_conv,'same'),bhv_matrix,'uni',false);
    hit_rate = cellfun(@(x) x(:,1)./(x(:,1)+x(:,4)),bhv_matrix_conv,'uni',false);
    false_alarm_rate = cellfun(@(x) x(:,3)./(x(:,3)+x(:,2)),bhv_matrix_conv,'uni',false);
    
    % Grab d'
    
    % get interpolated rate of conditions
    bhv_matrix_rate = cell(size(bhv_matrix));
    for curr_session = 1:length(bhv_matrix);
        for condition = 1:size(bhv_matrix{curr_session},2);
            
            curr_bin = bhv_matrix{curr_session}(:,condition);
            
            % if there's none of a current condition, rate is 0
            if ~any(curr_bin)
                bhv_matrix_rate{curr_session}(:,condition) = ...
                    zeros(length(curr_bin)-1,1);
                continue
            end
            
            curr_interp = cumsum(curr_bin);
            % set all zeros (except for the first and after last occurance) to NaN;
            curr_interp(find(curr_bin(2:find(curr_bin == 1,1,'last')) == 0) + 1) = NaN;
            curr_interp(isnan(curr_interp)) = ...
                interp1(find(~isnan(curr_interp)), ...
                curr_interp(~isnan(curr_interp)), ...
                find(isnan(curr_interp)));
            
            bhv_matrix_rate{curr_session}(:,condition) = diff(curr_interp);
        end
    end
    
    % smooth the rate (only in the forward direction for quick shifts)
    rate_smooth = 10;
    rate_smooth_kernel = ones(rate_smooth,1)/rate_smooth;
    bhv_matrix_rate_smooth_complete = cellfun(@(x) conv2(x,rate_smooth_kernel),bhv_matrix_rate,'uni',false);
    bhv_matrix_rate_smooth = cellfun(@(x) x(rate_smooth-1:end,:),bhv_matrix_rate_smooth_complete,'uni',false);
    
    % d' = zscore(hit rate (correct lick/correct lick + incorrect rejection)) -
    % zscore(false alarm rate (incorrect lick / incorrect lick + correct
    % rejection))
    bhv_raw_hit_falsealarm = cellfun(@(x) [x(:,1)./(x(:,1) + x(:,4)) ...
        x(:,3)./(x(:,3) + x(:,2))],bhv_matrix_rate_smooth,'uni',false);
    % d' can't use 0 or 1, adjust slightly
    bhv_raw_hit_falsealarm_adjust = cellfun(@(x) x - ((x == 1)*0.0001) + ...
        ((x == 0)*0.0001),bhv_raw_hit_falsealarm,'uni',false);
    bhv_dprime = cellfun(@(x) norminv(x(:,1),0,1) - norminv(x(:,2),0,1),bhv_raw_hit_falsealarm_adjust,'uni',false);
    
    % get chance lines by shuffling C/I lick/rejection trials
    num_shuff = 1000;
    bhv_dprime_shuff = cell(length(bhv_dprime),num_shuff);
    for curr_shuff = 1:num_shuff
        % shuffle C/I licks and C/I rejections
        bhv_matrix_unorder_shuffle = cellfun(@(curr_bhv) [shake(curr_bhv(:,[1 3]),2) ...
            shake(curr_bhv(:,[2 4]),2)],bhv_matrix,'uni',false);
        bhv_matrix_shuffle = cellfun(@(curr_bhv) curr_bhv(:,[1 3 2 4]), ...
            bhv_matrix_unorder_shuffle,'uni',false);

        % get interpolated rate of conditions
        bhv_matrix_rate = cell(size(bhv_matrix_shuffle));
        for curr_session = 1:length(bhv_matrix_shuffle);
            for condition = 1:size(bhv_matrix_shuffle{curr_session},2);
                
                curr_bin = bhv_matrix_shuffle{curr_session}(:,condition);
                
                % if there's none of a current condition, rate is 0
                if ~any(curr_bin)
                    bhv_matrix_rate{curr_session}(:,condition) = ...
                        zeros(length(curr_bin)-1,1);
                    continue
                end
                
                curr_interp = cumsum(curr_bin);
                % set all zeros (except for the first and after last occurance) to NaN;
                curr_interp(find(curr_bin(2:find(curr_bin == 1,1,'last')) == 0) + 1) = NaN;
                curr_interp(isnan(curr_interp)) = ...
                    interp1(find(~isnan(curr_interp)), ...
                    curr_interp(~isnan(curr_interp)), ...
                    find(isnan(curr_interp)));
                
                bhv_matrix_rate{curr_session}(:,condition) = diff(curr_interp);
            end           
        end
        
        % smooth the rate (only in the forward direction for quick shifts)
        bhv_matrix_rate_smooth_complete = cellfun(@(x) conv2(x,rate_smooth_kernel),bhv_matrix_rate,'uni',false);
        bhv_matrix_rate_smooth = cellfun(@(x) x(rate_smooth-1:end,:),bhv_matrix_rate_smooth_complete,'uni',false);
        
        % d' = zscore(hit rate (correct lick/correct lick + incorrect rejection)) -
        % zscore(false alarm rate (incorrect lick / incorrect lick + correct
        % rejection))
        bhv_raw_hit_falsealarm = cellfun(@(x) [x(:,1)./(x(:,1) + x(:,4)) ...
            x(:,3)./(x(:,3) + x(:,2))],bhv_matrix_rate_smooth,'uni',false);
        % d' can't use 0 or 1, adjust slightly
        bhv_raw_hit_falsealarm_adjust = cellfun(@(x) x - ((x == 1)*0.0001) + ...
            ((x == 0)*0.0001),bhv_raw_hit_falsealarm,'uni',false);
        bhv_dprime_shuff(:,curr_shuff) = cellfun(@(x) norminv(x(:,1),0,1) - norminv(x(:,2),0,1), ...
            bhv_raw_hit_falsealarm_adjust,'uni',false);
    end
    
    hit_rate_odor1cat = vertcat(hit_rate{:});
    hit_rate_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    hit_rate_odor2cat = vertcat(hit_rate{:});
    hit_rate_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    false_alarm_rate_odor1cat = vertcat(false_alarm_rate{:});
    false_alarm_rate_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    false_alarm_rate_odor2cat = vertcat(false_alarm_rate{:});
    false_alarm_rate_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    bhv_dprime_odor1cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    bhv_dprime_odor2cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    d_prime_chance = prctile(cell2mat(bhv_dprime_shuff),[2.5 97.5],2);
    
    figure;
    p1 = subplot(3,1,1); hold on;
    plot(hit_rate_odor1cat,'k');
    plot(hit_rate_odor2cat,'r');
    ylabel('Hit rate')
    title(['Behavior: ' animals{curr_animal}])
    xlim([0 sum(session_num_trials)]);
    ylim([-0.2 1.2])
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color',[0.5 0.5 0.5],'linestyle','--');
    end
    p2 = subplot(3,1,2); hold on;
    plot(false_alarm_rate_odor1cat,'k');
    plot(false_alarm_rate_odor2cat,'r');
    ylabel('False alarm rate')
    xlim([0 sum(session_num_trials)]);
    ylim([-0.2 1.2])
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color',[0.5 0.5 0.5],'linestyle','--');
    end
    p3 = subplot(3,1,3); hold on;
    plot(bhv_dprime_odor1cat,'k');
    plot(bhv_dprime_odor2cat,'r');
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color',[0.5 0.5 0.5],'linestyle','--');
    end
    plot(d_prime_chance,'b');
    xlim([0 sum(session_num_trials)]);
    ylabel('D''');
    xlabel('Trial')

end

%% Plot discrimination behavior (no interpolation, discrete)
% THIS IS WRONG!! see above cell

for curr_animal = 1;%:length(animals)
    
    rewarded_odor = cell(length(data_all(curr_animal).im),1);
    contingency_reversals = cell(length(data_all(curr_animal).im),1);
    bhv_matrix = cell(length(data_all(curr_animal).im),1);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        
        contingency_reversals{curr_session} = ...
            [0;diff(data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials))] ~= 0;
        
        bhv_matrix{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end
    
    session_num_trials = cellfun(@length,contingency_reversals);
    
    
    % Plot seperately: hit rate CL/(CL+IR), false alarm rate IL/(IL+CR)
    trial_range = 10; % surrounding trials to estimate by
    trial_conv = ones(trial_range,1);
    bhv_matrix_conv = cellfun(@(x) conv2(+x,trial_conv,'same'),bhv_matrix,'uni',false);
    hit_rate = cellfun(@(x) x(:,1)./(x(:,1)+x(:,4)),bhv_matrix_conv,'uni',false);
    false_alarm_rate = cellfun(@(x) x(:,3)./(x(:,3)+x(:,2)),bhv_matrix_conv,'uni',false);
    
    

    % d' = zscore(hit rate (correct lick/correct lick + incorrect rejection)) -
    % zscore(false alarm rate (incorrect lick / incorrect lick + correct
    % rejection))
    bhv_raw_hit_falsealarm = cellfun(@(x) [x(:,1)./(x(:,1) + x(:,4)) ...
        x(:,3)./(x(:,3) + x(:,2))],bhv_matrix_rate_smooth,'uni',false);
    bhv_all_mean = nanmean(vertcat(bhv_raw_hit_falsealarm{:}));
    bhv_all_std = nanstd(vertcat(bhv_raw_hit_falsealarm{:}));
    bhv_dprime = cellfun(@(x) diff(fliplr(bsxfun(@times,bsxfun(@minus,x, ...
        bhv_all_mean),bhv_all_std)),[],2),bhv_raw_hit_falsealarm,'uni',false);
    % this was to zscore by session
    %bhv_dprime = cellfun(@(x) diff(fliplr(bsxfun(@times,bsxfun(@minus,x, ...
    %    nanmean(x,1)),nanstd(x,[],1))),[],2),bhv_raw_hit_falsealarm,'uni',false);
    
    
    bhv_dprime_odor1cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor1cat(vertcat(rewarded_odor{:}) == 2) = NaN;
    bhv_dprime_odor2cat = vertcat(bhv_dprime{:});
    bhv_dprime_odor2cat(vertcat(rewarded_odor{:}) == 1) = NaN;
    
    figure; hold on;
    plot(bhv_dprime_odor1cat,'k');
    plot(bhv_dprime_odor2cat,'r');
    for i = cumsum(session_num_trials)
        line([i i],ylim,'color','k');
    end
    xlim([0 sum(session_num_trials)]);
    title(['Behavior: ' animals{curr_animal}])
    ylabel('D''');
    xlabel('Trial')
end




%% Plot licking behavior
correct_lick_trials_idx = cell(size(animals));
correct_rejection_trials_idx = cell(size(animals));
incorrect_lick_trials_idx = cell(size(animals));
incorrect_rejection_trials_idx = cell(size(animals));
lick_times = cell(size(animals));
answer_period_start = cell(size(animals));
reward_start = cell(size(animals));
for curr_animal = 1:length(animals)
    correct_lick_trials_idx{curr_animal} = cell(length(data_all(curr_animal).bhv),1);
    correct_rejection_trials_idx{curr_animal} = cell(length(data_all(curr_animal).bhv),1);
    incorrect_lick_trials_idx{curr_animal} = cell(length(data_all(curr_animal).bhv),1);
    incorrect_rejection_trials_idx{curr_animal} = cell(length(data_all(curr_animal).bhv),1);
    for curr_session = 1:length(data_all(curr_animal).bhv);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        correct_lick_trials_idx{curr_animal}{curr_session} = find(correct_lick_trials);
        correct_rejection_trials_idx{curr_animal}{curr_session} = find(correct_rejection_trials);
        incorrect_lick_trials_idx{curr_animal}{curr_session} = find(incorrect_lick_trials);
        incorrect_rejection_trials_idx{curr_animal}{curr_session} = find(incorrect_rejection_trials);
        
        lick_times{curr_animal}{curr_session} = cellfun(@(x) x.pokes.C(:,1) - x.states.apply_odor(1), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials),'uni',false);
        answer_period_start{curr_animal}{curr_session} = cellfun(@(x) x.states.answer(1) - x.states.apply_odor(1), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        discrimination_trials_idx = find(discrimination_trials);
        reward_start{curr_animal}{curr_session} = nan(sum(discrimination_trials),1);
        reward_start{curr_animal}{curr_session}(correct_lick_trials) = cellfun(@(x) x.states.correct_lick(1) - x.states.apply_odor(1), ...
            data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials_idx(correct_lick_trials)));
    end
end
    
figure; 
curr_animal = 3;
curr_session = 6;

subplot(2,2,1); hold on;
for i = 1:length(correct_lick_trials_idx{curr_animal}{curr_session})
   plot(lick_times{curr_animal}{curr_session}{correct_lick_trials_idx{curr_animal}{curr_session}(i)},i,'.k')
   plot(answer_period_start{curr_animal}{curr_session}(correct_lick_trials_idx{curr_animal}{curr_session}(i)),i,'.r');
   plot(reward_start{curr_animal}{curr_session}(correct_lick_trials_idx{curr_animal}{curr_session}(i)),i,'.g');
end
title('Correct lick')

subplot(2,2,3); hold on;
for i = 1:length(incorrect_lick_trials_idx{curr_animal}{curr_session})
   plot(lick_times{curr_animal}{curr_session}{incorrect_lick_trials_idx{curr_animal}{curr_session}(i)},i,'.k')
   plot(answer_period_start{curr_animal}{curr_session}(incorrect_lick_trials_idx{curr_animal}{curr_session}(i)),i,'.r');
end
title('Incorrect lick')

subplot(2,2,2); hold on;
for i = 1:length(correct_rejection_trials_idx{curr_animal}{curr_session})
   plot(answer_period_start{curr_animal}{curr_session}(correct_rejection_trials_idx{curr_animal}{curr_session}(i)),i,'.r');
   if isempty(lick_times{curr_animal}{curr_session}{correct_rejection_trials_idx{curr_animal}{curr_session}(i)})
       continue
   end
   plot(lick_times{curr_animal}{curr_session}{correct_rejection_trials_idx{curr_animal}{curr_session}(i)},i,'.k')
end
title('Correct rejection')

subplot(2,2,4); hold on;
for i = 1:length(incorrect_rejection_trials_idx{curr_animal}{curr_session})
   plot(answer_period_start{curr_animal}{curr_session}(incorrect_rejection_trials_idx{curr_animal}{curr_session}(i)),i,'.r');
   if isempty(lick_times{curr_animal}{curr_session}{incorrect_rejection_trials_idx{curr_animal}{curr_session}(i)})
       continue
   end
   plot(lick_times{curr_animal}{curr_session}{incorrect_rejection_trials_idx{curr_animal}{curr_session}(i)},i,'.k')
end
title('Incorrect rejection')


%% Lick analysis

% try to get times for separating 1) lick v rejection, 2) CL v IL
% get both: 
% - time (when the medians are significantly different) 
% - fraction correct classification (based on ideal observer)
% - how classifiable single trials are

%%%% CHANGE THIS! Make Conf intervals not overlap, not median/CI
%%% also the convolution doesn't take larger size into account, fix this


% Find each session when 1) CL/CR and 2) CL/IL diverge
lick_diverge_time = cell(size(animals));
correct_lick_diverge_time = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(data_all(curr_animal).im);
        % skip if no trials (happened when no discrim trials, rare)
        if isempty(lick_rate{curr_animal}{curr_session})
            continue
        end        
        
        % base divergence on median outside of other bootstrapped median
        % (only correct rejection trials, in case weird stuff incorrect)
        rejection_ci = bootci(1000,@nanmedian,lick_rate{curr_animal}{curr_session} ...
            ([correct_rejection_trials_idx{curr_animal}{curr_session}],:));
        % if small number of trials, just use range
        if length(incorrect_lick_trials_idx{curr_animal}{curr_session}) > 5
            incorrect_lick_ci = bootci(1000,@nanmedian,lick_rate{curr_animal}{curr_session} ...
                (incorrect_lick_trials_idx{curr_animal}{curr_session},:));
        else
            incorrect_lick_ci = [min(lick_rate{curr_animal}{curr_session} ...
                (incorrect_lick_trials_idx{curr_animal}{curr_session},:),[],1); ...
                max(lick_rate{curr_animal}{curr_session} ...
                (incorrect_lick_trials_idx{curr_animal}{curr_session},:),[],1)];
        end
        
        lick_median = nanmedian(lick_rate{curr_animal}{curr_session}( ...
            [correct_lick_trials_idx{curr_animal}{curr_session}; ...
            incorrect_lick_trials_idx{curr_animal}{curr_session}],:),1);
        
        correct_lick_median = nanmedian(lick_rate{curr_animal}{curr_session}( ...
            correct_lick_trials_idx{curr_animal}{curr_session},:),1);
        
        lick_diverge = correct_lick_median > rejection_ci(2,:);
        
        correct_lick_diverge = correct_lick_median < incorrect_lick_ci(1,:) | ...
            correct_lick_median > incorrect_lick_ci(2,:);
        
        % find first divergence point where different for certain time
        consecutive_divergence = 100; % in ms
        consecutive_divergence_samples = consecutive_divergence/1000*lick_resolution;
        
        lick_diverge_time{curr_animal}(curr_session) = (find(conv(+lick_diverge, ...
            ones(1,consecutive_divergence_samples),'valid')...
            == consecutive_divergence_samples,1))/lick_resolution;
        
        correct_lick_diverge_time{curr_animal}(curr_session) = (find(conv(+correct_lick_diverge, ...
            ones(1,consecutive_divergence_samples),'valid')...
            == consecutive_divergence_samples,1))/lick_resolution;
        
    end
    disp(curr_animal)
end

max_sessions = max(cellfun(@length,lick_diverge_time));
lick_diverge_time_pad = cellfun(@(x) padarray(x, ...
    [0 max_sessions-length(x)],nan,'post'),lick_diverge_time,'uni',false);
lick_diverge_time_cat = vertcat(lick_diverge_time_pad{:});

correct_lick_diverge_time_pad = cellfun(@(x) padarray(x, ...
    [0 max_sessions-length(x)],nan,'post'),correct_lick_diverge_time,'uni',false);
correct_lick_diverge_time_cat = vertcat(correct_lick_diverge_time_pad{:});

figure; 
subplot(2,1,1); hold on;
plot(lick_diverge_time_cat','color',[0.5 0.5 0.5]);
errorbar(nanmean(lick_diverge_time_cat),nanstd(lick_diverge_time_cat)./ ...
    sqrt(sum(~isnan(lick_diverge_time_cat))),'k','linewidth',2);
title('CL/CR diverge time')
subplot(2,1,2); hold on;
plot(correct_lick_diverge_time_cat','color',[0.5 0.5 0.5]);
errorbar(nanmean(correct_lick_diverge_time_cat),nanstd(correct_lick_diverge_time_cat)./ ...
    sqrt(sum(~isnan(correct_lick_diverge_time_cat))),'k','linewidth',2);
title('CL/IL diverge time')







%% Activity correlation

% split aligned activity into pyr/gad
aligned_df_gad = cellfun(@(x,y) cellfun(@(z) z(:,:,mice(y).interneurons),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);
aligned_df_pyr = cellfun(@(x,y) cellfun(@(z) z(:,:,setdiff(1:size(z,3),mice(y).interneurons)),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);


% get trial conditions and contingencies
condition_trials_all = cell(length(animals),1);
odor_trials_all = cell(length(animals),1);
for curr_animal = 1:length(animals)
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        condition_trials_all{curr_animal}{curr_session} = ...
            [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
        
        applied_odor_group = data_all(curr_animal).bhv(curr_session).applied_odor(discrimination_trials);
        rewarded_odor_group = data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        
        odor_trials_all{curr_animal}{curr_session} = [applied_odor_group rewarded_odor_group];
    end
end

% split activity further based on 1) applied odor, 2) trial condition
% possible applied odor = 1,2
% condition order: 1 = CL, 2 = CR, 3 = IL, 4 = IR

odor = [1 2];
condition = [1 2 3 4];

curr_activity = cellfun(@(x,y,z) cellfun(@(curr_df,curr_odor,curr_condition) curr_df( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),:,:), ...
    x,y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    aligned_df_pyr,odor_trials_all',condition_trials_all','uni',false);

curr_rewarded_odor = cellfun(@(y,z) cellfun(@(curr_odor,curr_condition) curr_odor( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),2), ...
    y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    odor_trials_all',condition_trials_all','uni',false);
    

curr_animal = 3;

% currently the mean of activity, ignoring timecourse
curr_activity_reshape = cellfun(@(x) ...
    reshape(permute(nanmean(x,2),[2 3 1]),[],size(x,1)),curr_activity{curr_animal},'uni',false);
cc_grid = cell(length(curr_activity_reshape));
for i = 1:length(curr_activity_reshape)
    for j = 1:length(curr_activity_reshape)
        curr_cc = corrcoef([curr_activity_reshape{i} curr_activity_reshape{j}]);
        cc_grid{i,j} = curr_cc(1:size(curr_activity_reshape{i},2),size(curr_activity_reshape{i},2)+1:end);
    end
end

cc_grid_cat = corrcoef(horzcat(curr_activity_reshape{:}));
session_trials_used = cumsum(cellfun(@(x) size(x,2),curr_activity_reshape));
figure; 
p1 = subplot(5,5,[7:10 12:15 17:20 22:25]);
imagesc(cc_grid_cat);colormap(gray);
curr_xlim = xlim;
p2 = subplot(5,5,2:5);
plot(vertcat(curr_rewarded_odor{curr_animal}{:}),'k');
xlim(curr_xlim);ylim([0 3]);axis off;
for i = 1:length(session_trials_used)
    line(repmat(session_trials_used(i),2,1),ylim,'color','r');
end
p3 = subplot(5,5,6:5:21);
plot(vertcat(curr_rewarded_odor{curr_animal}{:}),'k');
view(90,90);xlim(curr_xlim);ylim([0 3]);axis off;
for i = 1:length(session_trials_used)
    line(repmat(session_trials_used(i),2,1),ylim,'color','r');
end
linkaxes([p1 p2 p3],'x');


% compare running group of trials: extract diagonal strip
trials_compare = 100;
cc_grid_cat_compare = cc_grid_cat;
cc_grid_cat_compare(triu(true(size(cc_grid_cat))) | ...
    tril(true(size(cc_grid_cat)),-trials_compare)) = NaN;
figure; 
p1 = subplot(2,1,1);
plot(vertcat(curr_rewarded_odor{curr_animal}{:}),'k','linewidth',2);
ylim([0 3]);
for i = 1:length(session_trials_used)
    line(repmat(session_trials_used(i),2,1),ylim,'color','r');
end
ylabel('Rewarded odor')
p2 = subplot(2,1,2);
plot(smooth(nanmean(cc_grid_cat_compare),10),'k','linewidth',2);
ylabel('Activity correlation')
linkaxes([p1 p2],'x');


%%%%% Same as before, but comparing two groups instead of within groups
% condition order: 1 = CL, 2 = CR, 3 = IL, 4 = IR

odor = [1];
condition = [1];

curr_activity1 = cellfun(@(x,y,z) cellfun(@(curr_df,curr_odor,curr_condition) curr_df( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),:,:), ...
    x,y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    aligned_df_pyr,odor_trials_all',condition_trials_all','uni',false);

curr_rewarded_odor1 = cellfun(@(y,z) cellfun(@(curr_odor,curr_condition) curr_odor( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),2), ...
    y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    odor_trials_all',condition_trials_all','uni',false);
    
odor = [2];
condition = [1];

curr_activity2 = cellfun(@(x,y,z) cellfun(@(curr_df,curr_odor,curr_condition) curr_df( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),:,:), ...
    x,y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    aligned_df_pyr,odor_trials_all',condition_trials_all','uni',false);

curr_rewarded_odor2 = cellfun(@(y,z) cellfun(@(curr_odor,curr_condition) curr_odor( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),2), ...
    y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    odor_trials_all',condition_trials_all','uni',false);

curr_animal = 3;

% currently the mean of activity, ignoring timecourse
curr_activity_reshape1 = cellfun(@(x) ...
    reshape(permute(nanmean(x,2),[2 3 1]),[],size(x,1)),curr_activity1{curr_animal},'uni',false);

curr_activity_reshape2 = cellfun(@(x) ...
    reshape(permute(nanmean(x,2),[2 3 1]),[],size(x,1)),curr_activity2{curr_animal},'uni',false);

cc_grid = cell(length(curr_activity_reshape1),length(curr_activity_reshape2));
for i = 1:length(curr_activity_reshape1)
    for j = 1:length(curr_activity_reshape2)
        curr_cc = corrcoef([curr_activity_reshape1{i} curr_activity_reshape2{j}]);
        cc_grid{i,j} = curr_cc(1:size(curr_activity_reshape1{i},2),size(curr_activity_reshape1{i},2)+1:end);
    end
end

cc_grid_cat_full = corrcoef([horzcat(curr_activity_reshape1{:}) ...
    horzcat(curr_activity_reshape2{:})]);
cc_grid_cat = cc_grid_cat_full(1:size(horzcat(curr_activity_reshape1{:}),2), ...
    size(horzcat(curr_activity_reshape1{:}),2)+1:end);
session_trials_used1 = cumsum(cellfun(@(x) size(x,2),curr_activity_reshape1));
session_trials_used2 = cumsum(cellfun(@(x) size(x,2),curr_activity_reshape2));
figure; 
p1 = subplot(5,5,[7:10 12:15 17:20 22:25]);
imagesc(cc_grid_cat);colormap(gray);
curr_lim1 = ylim;
curr_lim2 = xlim;
p3 = subplot(5,5,6:5:21);
plot(vertcat(curr_rewarded_odor1{curr_animal}{:}),'k');
view(90,90);xlim(curr_lim1);ylim([0 3]);axis off;
for i = 1:length(session_trials_used1)
    line(repmat(session_trials_used1(i),2,1),ylim,'color','r');
end
p2 = subplot(5,5,2:5);
plot(vertcat(curr_rewarded_odor2{curr_animal}{:}),'k');
xlim(curr_lim2);ylim([0 3]);axis off;
for i = 1:length(session_trials_used2)
    line(repmat(session_trials_used2(i),2,1),ylim,'color','r');
end





%%% do the last thing, but combine animals
odor = [1];
condition = [1];

curr_activity1 = cellfun(@(x,y,z) cellfun(@(curr_df,curr_odor,curr_condition) curr_df( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),:,:), ...
    x,y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    aligned_df_pyr,odor_trials_all',condition_trials_all','uni',false);

curr_rewarded_odor1 = cellfun(@(y,z) cellfun(@(curr_odor,curr_condition) curr_odor( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),2), ...
    y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    odor_trials_all',condition_trials_all','uni',false);
    
odor = [2];
condition = [1];

curr_activity2 = cellfun(@(x,y,z) cellfun(@(curr_df,curr_odor,curr_condition) curr_df( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),:,:), ...
    x,y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    aligned_df_pyr,odor_trials_all',condition_trials_all','uni',false);

curr_rewarded_odor2 = cellfun(@(y,z) cellfun(@(curr_odor,curr_condition) curr_odor( ...
    ismember(curr_odor(:,1),odor) & any(curr_condition(:,condition),2),2), ...
    y(cellfun(@(r) ~isempty(r),y)), ...
    z(cellfun(@(r) ~isempty(r),z)),'uni',false), ...
    odor_trials_all',condition_trials_all','uni',false);



max_sessions = max(cellfun(@length,curr_activity1));
cc_grid = cell(max_sessions,max_sessions,length(animals));

for curr_animal = 1:length(animals)
    
    % currently the mean of activity, ignoring timecourse
    curr_activity_reshape1 = cellfun(@(x) ...
        reshape(permute(nanmean(x,2),[2 3 1]),[],size(x,1)),curr_activity1{curr_animal},'uni',false);
    
    curr_activity_reshape2 = cellfun(@(x) ...
        reshape(permute(nanmean(x,2),[2 3 1]),[],size(x,1)),curr_activity2{curr_animal},'uni',false);
    
    for i = 1:length(curr_activity_reshape1)
        for j = 1:length(curr_activity_reshape2)
            curr_cc = corrcoef([curr_activity_reshape1{i} curr_activity_reshape2{j}]);
            cc_grid{i,j,curr_animal} =  ...
                curr_cc(1:size(curr_activity_reshape1{i},2),size(curr_activity_reshape1{i},2)+1:end);
        end
    end
    
end


max_sessions = max(cellfun(@length,curr_activity1));
cc_grid = cell(max_sessions,max_sessions,length(animals));

for curr_animal = 1:length(animals)
    
    % currently the mean timecourse of activity
    curr_activity_reshape1 = cellfun(@(x) ...
        reshape(permute(x,[2 3 1]),[],size(x,1)),curr_activity1{curr_animal},'uni',false);
    
    curr_activity_reshape2 = cellfun(@(x) ...
        reshape(permute(x,[2 3 1]),[],size(x,1)),curr_activity2{curr_animal},'uni',false);
    
    for i = 1:length(curr_activity_reshape1)
        for j = 1:length(curr_activity_reshape2)
            curr_cc = corrcoef([curr_activity_reshape1{i} curr_activity_reshape2{j}]);
            cc_grid{i,j,curr_animal} =  ...
                curr_cc(1:size(curr_activity_reshape1{i},2),size(curr_activity_reshape1{i},2)+1:end);
        end
    end
    
end


%% Get cells significantly modulated during lick trials (CL+IL vs R)
% the code cell below this is better statistics but longer, this is quick

lick_p = cell(size(animals));
for curr_animal = 1:length(animals)
    lick_p{curr_animal} = cell(length(data_all(curr_animal).im),1);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        % get trials with 'discrimination' paradigm that were also imaged for the
        % course of the epoch time
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
         
         % get number of usable trials and cells
         num_trials = sum(discrimination_trials);
         num_cells = size(data_all(curr_animal).im(curr_session).roi_trace_df{1},1);
         
         if num_trials == 0
             lick_p{curr_animal}{curr_session} = nan(1,num_cells);
             continue
         end

        % take aligned time from above
        
        
        % classify significantly different cells based on shuffled F-statistic of
        % mean df during trials.
        odor_aligned_df_mean = permute(nanmean(aligned_df_all{curr_animal}{curr_session},2),[1 3 2]);

        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        
        % NOTE! this is not the correct way to do these statistics, but
        % they're probably close and good enough for preliminary analysis
        lick_trials = correct_lick_trials | incorrect_lick_trials;
        lick_p{curr_animal}{curr_session} = arrayfun(@(x) ranksum(odor_aligned_df_mean(correct_lick_trials,x), ...
            odor_aligned_df_mean(correct_rejection_trials,x)),1:size(odor_aligned_df_mean,2));

    end
    disp(curr_animal)
end

% note!! changed condition_sig_cells to p value, different name
lick_sig_cells = cellfun(@(x) cellfun(@(y) y < 0.025,x,'uni',false), ...
    lick_p,'uni',false);

lick_sig_gad = cellfun(@(x,mouse_idx) cellfun(@(y) ...
    y(mice(mouse_idx).interneurons), ...
    x,'uni',false), ...
    lick_sig_cells,mice_idx,'uni',false);

lick_sig_pyr = cellfun(@(x,mouse_idx) cellfun(@(y) ...
    y(setdiff(1:length(y),mice(mouse_idx).interneurons)), ...
    x,'uni',false), ...
    lick_sig_cells,mice_idx,'uni',false);

% plot lick population overlap
max_sessions = max(cellfun(@length,lick_sig_pyr));
lick_pop_overlap = nan(max_sessions,max_sessions,length(animals));
for curr_animal = 1:length(animals)
   curr_lick_cells =  vertcat(lick_sig_pyr{curr_animal}{:});
   curr_lick_overlap = +curr_lick_cells*+curr_lick_cells';
   lick_pop_overlap(1:size(curr_lick_cells,1), ...
       1:size(curr_lick_cells,1),curr_animal) = curr_lick_overlap./ ...
       sqrt(repmat(sum(curr_lick_cells,2),1,size(curr_lick_cells,1)).* ...
       repmat(sum(curr_lick_cells,2),1,size(curr_lick_cells,1))');
end

% fraction lick cells over time
frac_lick_gad = cellfun(@(x) cellfun(@nanmean,x),lick_sig_gad,'uni',false);
frac_lick_pyr = cellfun(@(x) cellfun(@nanmean,x),lick_sig_pyr,'uni',false);


%% Get divergence time for cells which are CR-modulated related



activity_diverge_time_pyr = cell(size(animals));
for curr_animal = 1:length(animals)
    activity_diverge_time_pyr{curr_animal} = cell(length(lick_sig_cells{curr_animal}),1);
    for curr_session = 1:length(lick_sig_cells{curr_animal});
        
        curr_gad = false(size(data_all(curr_animal).im(curr_session).roi_trace_df{1},1),1);
        curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
        curr_pyr = ~curr_gad;
        
        curr_cells = find(lick_sig_cells{curr_animal}{curr_session} & curr_pyr');
        
        activity_diverge_time_pyr{curr_animal}{curr_session} = ...
            nan(length(curr_pyr),1);
        
        CL_ci = arrayfun(@(x) bootci(100,@nanmean,aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1),:,x)), ...
            curr_cells,'uni',false);
        
        CR_ci = arrayfun(@(x) bootci(100,@nanmean,aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,2),:,x)), ...
            curr_cells,'uni',false);
        
        activity_diverge = cellfun(@(x,y) x(1,:) > y(2,:) | ...
            x(2,:) < y(1,:),CL_ci,CR_ci,'uni',false);
        
        % find first divergence point where different for certain time
        consecutive_divergence = 3; % in frames
        
        curr_divergence = cellfun(@(x) ... 
            find(conv(+x,ones(1,consecutive_divergence),'valid') == ...
            consecutive_divergence,1),activity_diverge,'uni',false); 
        
        use_curr_divergence = cellfun(@(x) ~isempty(x),curr_divergence);
        
        activity_diverge_time_pyr{curr_animal}{curr_session}( ...
            curr_cells(use_curr_divergence)) = ...
            vertcat(curr_divergence{use_curr_divergence});
    end
    disp(curr_animal)
end








%% Get significantly modulated cells
% Note: animals are allowed to lick during odor time
condition_p = cell(size(animals));
for curr_animal = 1:length(animals)
    condition_p{curr_animal} = cell(length(data_all(curr_animal).im),1);
    for curr_session = 1:length(data_all(curr_animal).im);
        
%         % user settings
%         epoch_time = [0 2]; % time to analyze surrounding event (sec)
%         loop_frames = 1000;
         num_epochs = 1;
%         
%         %%%%%%%%%%%%%%%%
%         
%         epoch_frames = round(epoch_time * data_all(curr_animal).im(curr_session).framerate);
%         
        % get trials with 'discrimination' paradigm that were also imaged for the
        % course of the epoch time
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
         
         % get number of usable trials and cells
         num_trials = sum(discrimination_trials);
         num_cells = size(data_all(curr_animal).im(curr_session).roi_trace_df{1},1);
         
         if num_trials == 0
             condition_p{curr_animal}{curr_session} = nan(5,num_cells);
             continue
         end
%         
%         % skip if number of trials is 0
%         if num_trials == 0
%             condition_sig_cells{curr_animal}{curr_session} = nan(5,num_cells);
%             continue
%         end
%         
%         % pull out the loop in which usable trials occurred
%         discrimination_trials_loop = data_all(curr_animal).bhv(curr_session).xsg_trials(discrimination_trials);
%         
%         % get the start time of odor onset and answer period (to nearest frame)
%         odor_onsets = cellfun(@(x) round(x.states.apply_odor(1)), ...
%             data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
%         answer_onsets = cellfun(@(x) round(x.states.answer(1)), ...
%             data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
%         
%         % get the imaging data surrounding odor/answer epochs
%         odor_aligned_df = permute(cell2mat(arrayfun(@(curr_trial) ...
%             data_all(curr_animal).im(curr_session).roi_trace_df{discrimination_trials_loop(curr_trial)} ...
%             (:,odor_onsets(curr_trial) + epoch_frames(1) : ...
%             odor_onsets(curr_trial) + epoch_frames(2)), permute(1:num_trials,[1 3 2]),'uni',false)),[3 2 1]);
%         
%         answer_aligned_df = permute(cell2mat(arrayfun(@(curr_trial) ...
%             data_all(curr_animal).im(curr_session).roi_trace_df{discrimination_trials_loop(curr_trial)} ...
%             (:,answer_onsets(curr_trial) + epoch_frames(1) : ...
%             answer_onsets(curr_trial) + epoch_frames(2)), permute(1:num_trials,[1 3 2]),'uni',false)),[3 2 1]);
                
        % take aligned time from above
        
        
        % classify significantly different cells based on shuffled F-statistic of
        % mean df during trials.
        odor_aligned_df_mean = permute(nanmean(aligned_df_all{curr_animal}{curr_session},2),[1 3 2]);

        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        
        
        % define anova groups: 1) epoch, 2) odor, 3) lick, 4) rewarded
        
        epoch_group = cell2mat(arrayfun(@(x) x*ones(num_trials,1),1:num_epochs,'uni',false));
        applied_odor_group = data_all(curr_animal).bhv(curr_session).applied_odor(discrimination_trials);
        rewarded_odor_group = data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);
        lick_group = correct_lick_trials | incorrect_lick_trials;
        rewarded_group = correct_lick_trials;
        
        anova_groups = {epoch_group,...
            repmat(applied_odor_group,num_epochs,1), ...
            repmat(rewarded_odor_group,num_epochs,1), ...
            repmat(lick_group,num_epochs,1), ...
            repmat(rewarded_group,num_epochs,1)};
        
        % run straight
        %[p t] = arrayfun(@(curr_cell) ...
        %    anovan([odor_aligned_df_mean(:,curr_cell);answer_aligned_df_mean(:,curr_cell)], ...
        %    anova_groups,'display','off'),1:num_cells,'uni',false);
        % Using entire odor onset to (would-be) answer offset
        [p t] = arrayfun(@(curr_cell) ...
            anovan([odor_aligned_df_mean(:,curr_cell)], ...
            anova_groups,'display','off'),1:num_cells,'uni',false);
        
        f_stat = cellfun(@(x) vertcat(x{2:end,6}),t,'uni',false);
        
        % run shuffle anova
        num_shuffle = 100;
        f_stat_shuffle = cell(num_cells,num_shuffle);
        for curr_shuff = 1:num_shuffle
            curr_anova_groups_shuffle = cellfun(@(x) shake(x),anova_groups,'uni',false);
            
            %[p t] = arrayfun(@(curr_cell) ...
            %    anovan([odor_aligned_df_mean(:,curr_cell);answer_aligned_df_mean(:,curr_cell)], ...
            %    curr_anova_groups_shuffle,'display','off'),1:num_cells,'uni',false);
            % Using entire odor onset to (would-be) answer offset
            [p t] = arrayfun(@(curr_cell) ...
                anovan([odor_aligned_df_mean(:,curr_cell)], ...
                curr_anova_groups_shuffle,'display','off'),1:num_cells,'uni',false);
            
            f_stat_shuffle(:,curr_shuff) = cellfun(@(x) vertcat(x{2:end,6}),t,'uni',false);
        end
        
        f_stat_shuffle_rank = arrayfun(@(curr_cell) ...
            tiedrank([horzcat(f_stat_shuffle{curr_cell,:}),f_stat{curr_cell}]')',1:num_cells,'uni',false);
        
        f_stat_p = cellfun(@(x) x(:,end)./size(x,2),f_stat_shuffle_rank,'uni',false);
        
        condition_p{curr_animal}{curr_session} = horzcat(f_stat_p{:});
        disp(curr_session);
    end
end

% note!! changed condition_sig_cells to p value, different name
condition_sig_cells = cellfun(@(x) cellfun(@(y) y > 0.975 | y < 0.025,x,'uni',false), ...
    condition_p,'uni',false);

% fraction lick cells over time
frac_lick_gad = cellfun(@(x,y) cellfun(@(z) nanmean(z(4,mice(y).interneurons), ...
    2),x),condition_sig_cells,mice_idx,'uni',false);
frac_lick_pyr = cellfun(@(x,y) cellfun(@(z) nanmean(z(4,setdiff(1:size(z,2),mice(y).interneurons)), ...
    2),x),condition_sig_cells,mice_idx,'uni',false);

lick_pyr = cellfun(@(x,y) cellfun(@(z) z(4,setdiff(1:size(z,2),mice(y).interneurons)), ...
    x,'uni',false),condition_sig_cells,mice_idx,'uni',false);
lick_popcorr = cellfun(@(x) corrcoef(vertcat(x{:})'),lick_pyr,'uni',false);

frac_contingency_pyr = cellfun(@(x,y) cellfun(@(z) nanmean(z(3,setdiff(1:size(z,2),mice(y).interneurons)), ...
    2),x),condition_sig_cells,mice_idx,'uni',false);



for i = 1:6
   figure; hold on
   plot(frac_lick_pyr{i},'k');
   plot(frac_lick_gad{i},'r');  
   title(animals{i})
end



%%  testing ground: trial-trial PCA

% PCA on trials: does O2CL approach original O1CL, then both diverge?
rewarded_odor = cell(size(animals));
bhv_matrix = cell(size(animals));
for curr_animal = 1:length(animals);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_animal}{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);

        bhv_matrix{curr_animal}{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end    
end


CL1_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL1_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1) ...
            & rewarded_odor{curr_animal}{curr_session} == 1,:,:),1),[3 2 1])';
    end
end

CL2_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL2_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1) ...
            & rewarded_odor{curr_animal}{curr_session} == 2,:,:),1),[3 2 1])';
    end
end

curr_animal = 2;

curr_gad = false(size(CL1_aligned_df_mean{curr_animal}{1},2),1);
curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
curr_pyr = ~curr_gad;

CL1_use_sessions = cellfun(@(x) ~all(isnan(x(:))),CL1_aligned_df_mean{curr_animal});
CL2_use_sessions = cellfun(@(x) ~all(isnan(x(:))),CL2_aligned_df_mean{curr_animal});

cat_df = [vertcat(CL1_aligned_df_mean{curr_animal}{CL1_use_sessions}); ...
    vertcat(CL2_aligned_df_mean{curr_animal}{CL2_use_sessions})];

cat_df_pyr = cat_df(:,curr_pyr);

% normalize by z-scoring
cat_df_pyr_norm = zscore(cat_df_pyr);
[coeff score latent] = princomp(cat_df_pyr_norm);

% split scores back up into contingency/day
CL1_scores = score(1:size(vertcat(CL1_aligned_df_mean{curr_animal}{CL1_use_sessions}),1),1:3);
CL2_scores = score(size(vertcat(CL1_aligned_df_mean{curr_animal}{CL1_use_sessions}),1)+1:end,1:3);

CL1_scores_sessions = mat2cell(CL1_scores,length(epoch_frames(1):epoch_frames(2))* ...
    ones(sum(CL1_use_sessions),1),3);
CL2_scores_sessions = mat2cell(CL2_scores,length(epoch_frames(1):epoch_frames(2))* ...
    ones(sum(CL2_use_sessions),1),3);

figure; hold on;
col = hot(15);
for i = 1:length(CL1_scores_sessions)
    plot3(CL1_scores_sessions{i}(:,1),CL1_scores_sessions{i}(:,2), ...
        CL1_scores_sessions{i}(:,3),'color',col(i,:));    
end
col = hot(15);
for i = 1:length(CL2_scores_sessions)
    plot3(CL2_scores_sessions{i}(:,1),CL2_scores_sessions{i}(:,2), ...
        CL2_scores_sessions{i}(:,3),'color',col(i,:),'linestyle','--');    
end




% 
% % apply coeffs to data (why is this different from scores?)
% CL1_aligned_df_mean_pca = cellfun(@(df) df(:,curr_pyr)*coeff(:,1:3), ...
%     CL1_aligned_df_mean{curr_animal},'uni',false);
% CL2_aligned_df_mean_pca = cellfun(@(df) df(:,curr_pyr)*coeff(:,1:3), ...
%     CL2_aligned_df_mean{curr_animal},'uni',false);
% 







%% Testing: PCA variance over time

% PCA on trials: does O2CL approach original O1CL, then both diverge?
rewarded_odor = cell(size(animals));
bhv_matrix = cell(size(animals));
for curr_animal = 1:length(animals);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_animal}{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);

        bhv_matrix{curr_animal}{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end    
end

% pull out CL trials
CL_aligned_df = cellfun(@(x,y) cellfun(@(a,b) a(b(:,1),:,:),x,y, ...
    'uni',false),aligned_df_all,bhv_matrix,'uni',false);
CL_aligned_df_reshape = cellfun(@(x) cellfun(@(y) ...
    reshape(permute(y,[2 1 3]), ...
    [],size(y,3)),x,'uni',false),CL_aligned_df,'uni',false);

CL_aligned_df_trial_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL_aligned_df_trial_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1),:,:),2),[3 1 2])';
    end
end

CL_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1),:,:),1),[3 2 1])';
    end
end

% pull out CR trials
CR_aligned_df = cellfun(@(x,y) cellfun(@(a,b) a(b(:,2),:,:),x,y, ...
    'uni',false),aligned_df_all,bhv_matrix,'uni',false);
CR_aligned_df_reshape = cellfun(@(x) cellfun(@(y) ...
    reshape(permute(y,[2 1 3]), ...
    [],size(y,3)),x,'uni',false),CR_aligned_df,'uni',false);

CR_aligned_df_trial_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CR_aligned_df_trial_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,2),:,:),2),[3 1 2])';
    end
end

CR_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CR_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,2),:,:),1),[3 2 1])';
    end
end


for curr_animal = 1:length(animals)
    
    curr_gad = false(size(CL1_aligned_df_mean{curr_animal}{end},2),1);
    curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
    curr_pyr = ~curr_gad;
    
    curr_data = [CL_aligned_df_mean{curr_animal};CR_aligned_df_mean{curr_animal}];
    
    % normalize by z-scoring
    pyr_mean_norm = cellfun(@(x,y) zscore([x(:,curr_pyr);y(:,curr_pyr)]),...
        CL_aligned_df_mean{curr_animal},CR_aligned_df_mean{curr_animal},'uni',false);

    % get the PCs
    [coeff score latent] = cellfun(@princomp,pyr_mean_norm,'uni',false);
    pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);
end






%% EVERYTHING BELOW: EXPERIMENTAL PCA STUFF, REQUIRES THIS CELL

% PCA on trials: does O2CL approach original O1CL, then both diverge?
rewarded_odor = cell(size(animals));
bhv_matrix = cell(size(animals));
for curr_animal = 1:length(animals);
    for curr_session = 1:length(data_all(curr_animal).im);
        
        discrimination_trials = cellfun(@(x,y) ...
            strcmp('Discrimination',x) & ...
            y.states.bitcode(1) > 1 & ...
            y.states.iti(2) < loop_frames, ...
            data_all(curr_animal).bhv(curr_session).session_type,data_all(curr_animal).bhv(curr_session).bhv_frames) &...
            data_all(curr_animal).bhv(curr_session).imaged_trials;
        
        % get trial conditions
        correct_lick_trials = cellfun(@(x) isfield(x.states,'correct_lick') && ...
            ~isempty(x.states.correct_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        correct_rejection_trials = cellfun(@(x) isfield(x.states,'correct_rejection') && ...
            ~isempty(x.states.correct_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_lick_trials = cellfun(@(x) isfield(x.states,'incorrect_lick') && ...
            ~isempty(x.states.incorrect_lick),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        incorrect_rejection_trials = cellfun(@(x) isfield(x.states,'incorrect_rejection') && ...
            ~isempty(x.states.incorrect_rejection),data_all(curr_animal).bhv(curr_session).bhv_frames(discrimination_trials));
        
        rewarded_odor{curr_animal}{curr_session} = ...
            data_all(curr_animal).bhv(curr_session).rewarded_odor(discrimination_trials);

        bhv_matrix{curr_animal}{curr_session} = [correct_lick_trials correct_rejection_trials ...
            incorrect_lick_trials incorrect_rejection_trials];
    end    
end

CL_aligned_df = cellfun(@(x,y) cellfun(@(a,b) a(b(:,1),:,:),x,y, ...
    'uni',false),aligned_df_all,bhv_matrix,'uni',false);
CL_aligned_df_reshape = cellfun(@(x) cellfun(@(y) ...
    reshape(permute(y,[2 1 3]), ...
    [],size(y,3)),x,'uni',false),CL_aligned_df,'uni',false);

CL_aligned_df_trial_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL_aligned_df_trial_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1),:,:),2),[3 1 2])';
    end
end

CL_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CL_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,1),:,:),1),[3 2 1])';
    end
end

CR_aligned_df_mean = cell(size(animals));
for curr_animal = 1:length(animals)
    for curr_session = 1:length(aligned_df_all{curr_animal});
        CR_aligned_df_mean{curr_animal}{curr_session} = ...
            permute(nanmean(aligned_df_all{curr_animal}{curr_session}(...
            bhv_matrix{curr_animal}{curr_session}(:,2),:,:),1),[3 2 1])';
    end
end

%%

curr_animal = 6;

curr_gad = false(size(CL_aligned_df_mean{curr_animal}{end},2),1);
curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
curr_pyr = ~curr_gad;


% combine L/R normalize by z-scoring
CL_CR_pyr = cellfun(@(x,y,z) [x(:,curr_pyr & z');y(:,curr_pyr & z')], ...
    CL_aligned_df_mean{curr_animal}, ...
    CR_aligned_df_mean{curr_animal}, ...
    lick_sig_cells{curr_animal}','uni',false);
CL_CR_pyr_norm = cellfun(@(x) zscore(x),CL_CR_pyr,'uni',false);

[coeff score latent] = cellfun(@princomp,CL_CR_pyr_norm,'uni',false);

pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);



%% 

% PROB INCLUDE LICKING HERE TOO? The variance explained of licking should always
% be the same
pc_mean_var_cat_all = cell(size(animals));
pc_trial_var_cat_all = cell(size(animals));
pc_trial_mean_var_cat_all = cell(size(animals));

figure;
for curr_animal = 1:length(animals)

curr_gad = false(size(CL1_aligned_df_mean{curr_animal}{end},2),1);
curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
curr_pyr = ~curr_gad;

% normalize by z-scoring
CL_pyr_trial_norm = cellfun(@(x) zscore(x(:,curr_pyr)),CL_aligned_df_reshape{curr_animal},'uni',false);
CL_pyr_trial_mean_norm = cellfun(@(x) zscore(x(:,curr_pyr)),CL_aligned_df_trial_mean{curr_animal},'uni',false);
CL_pyr_norm = cellfun(@(x) zscore(x(:,curr_pyr)),CL_aligned_df_mean{curr_animal},'uni',false);


% get the variance/PCs of mean
CL_pyr_norm_var = cellfun(@(x) nanmean(var(x,[],2)),CL_pyr_norm);

[coeff score latent] = cellfun(@princomp,CL_pyr_norm,'uni',false);
pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);

pc_mean_var_cat = horzcat(pc_variance_explained{:});
pc_mean_var_cat_all{curr_animal} = pc_mean_var_cat;
subplot(length(animals),3,curr_animal*3-2);hold on
plot(nanmean(pc_mean_var_cat(:,1:3),2),'k')
plot(nanmean(pc_mean_var_cat(:,4:6),2),'r')
plot(nanmean(pc_mean_var_cat(:,7:10),2),'b')
xlabel('Number of principle components')
ylabel('Cumulative explained variance')
title('Cells (mean)')


% get variance/PCs of trial-trial
CL_pyr_trial_norm_var = cellfun(@(x) nanmean(var(x,[],2)),CL_pyr_trial_norm);

[coeff score latent] = cellfun(@princomp,CL_pyr_trial_norm,'uni',false);
pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);

pc_trial_var_cat = horzcat(pc_variance_explained{:});
pc_trial_var_cat_all{curr_animal} = pc_trial_var_cat;
subplot(length(animals),3,curr_animal*3-1);hold on
plot(nanmean(pc_trial_var_cat(:,1:3),2),'k')
plot(nanmean(pc_trial_var_cat(:,4:6),2),'r')
plot(nanmean(pc_trial_var_cat(:,7:10),2),'b')
xlabel('Number of principle components')
ylabel('Cumulative explained variance')
title('Cells (trial-trial)')


% get variance/PCs of mean trial-trial (for groups of cross-category cells)
CL_pyr_trial_mean_norm_var = cellfun(@(x) nanmean(var(x,[],2)),CL_pyr_trial_mean_norm);

[coeff score latent] = cellfun(@princomp,CL_pyr_trial_mean_norm,'uni',false);
pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);

pc_trial_mean_var_cat = horzcat(pc_variance_explained{:});
pc_trial_mean_var_cat_all{curr_animal} = pc_trial_mean_var_cat;
subplot(length(animals),3,curr_animal*3);hold on
plot(nanmean(pc_trial_mean_var_cat(:,1:3),2),'k')
plot(nanmean(pc_trial_mean_var_cat(:,4:6),2),'r')
plot(nanmean(pc_trial_mean_var_cat(:,7:10),2),'b')
xlabel('Number of principle components')
ylabel('Cumulative explained variance')
title('Cells (mean trial-trial)')

end

% get correlation within PCs by PC contribution
pc_bins = [-1:0.2:0.8 Inf];
pc_corr = cell(length(coeff),1);
for curr_session = 1:length(coeff)
    pc_corr{curr_session} = nan(size(coeff{curr_session},2),length(pc_bins));
   for curr_pc = 1:size(coeff{curr_session},2)
       % get the trial-trial correlation of normalized activity
       pyr_coeff = corrcoef(CL_pyr_trial_norm{curr_session});
       
       % normalize the pcs (1 to -1)
       pc = coeff{curr_session}(:,curr_pc);
       norm_pc = pc./max(abs(pc));
       
       % bin cells by pc contribution, pull out correlation
       [n bins] = histc(norm_pc,pc_bins);
       for curr_bin = 1:length(pc_bins)
           pc_corr{curr_session}(curr_pc,curr_bin) = ...
               nanmedian(AP_itril(pyr_coeff(bins == curr_bin, ...
               bins == curr_bin),-1));
       end
   end
end




%CL1_use_sessions = cellfun(@(x) ~all(isnan(x(:))),CL1_aligned_df_mean{curr_animal});
%CL2_use_sessions = cellfun(@(x) ~all(isnan(x(:))),CL2_aligned_df_mean{curr_animal});





%%

curr_animal = 2;
curr_session = 10;
curr_gad = false(size(CL1_aligned_df_mean{curr_animal}{end},2),1);
curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
curr_pyr = ~curr_gad;

CL_pyr_trial_mean_norm = zscore(CL_aligned_df_trial_mean{curr_animal}{curr_session}(:,curr_pyr));

[coeff score latent] = princomp(CL_pyr_trial_mean_norm);
pc_variance_explained = cumsum(latent)./sum(latent);

num_shuff = 1000;
pc_variance_explained_shuff = nan(length(pc_variance_explained),num_shuff);
for i = 1:num_shuff
    [coeff score latent] = princomp(shake(CL_pyr_trial_mean_norm,1));
    pc_variance_explained_shuff(:,i) = cumsum(latent)./sum(latent);
    disp(i);
end



%%

curr_animal = 1;
curr_gad = false(size(CL1_aligned_df_mean{curr_animal}{end},2),1);
curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
curr_pyr = ~curr_gad;

CL_pyr_trial_mean_norm = cellfun(@(x) zscore(x(:,curr_pyr)),CL_aligned_df_trial_mean{curr_animal},'uni',false);

[coeff score latent] = cellfun(@princomp,CL_pyr_trial_mean_norm,'uni',false);
pc_variance_explained = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);


num_shuff = 1000;
pc_variance_explained_shuff = cell(num_shuff,length(pc_variance_explained));
for i = 1:num_shuff
    [coeff score latent] = cellfun(@(x) princomp(shake(x,1)),CL_pyr_trial_mean_norm,'uni',false);
    pc_variance_explained_shuff(i,:) = cellfun(@(x) cumsum(x)./sum(x),latent,'uni',false);
    disp(i);
end


pc_variance_explained_ci = cell(1,length(pc_variance_explained));
for i = 1:length(pc_variance_explained)
    pc_variance_explained_ci{i} = prctile(horzcat(pc_variance_explained_shuff{:,i}),[2.5 97.5],2);
end



figure;
for i = 1:10
   subplot(2,5,i); hold on
   plot(pc_variance_explained_ci{i},'r');
   plot(pc_variance_explained{i},'k','linewidth',2);
end



%% Get fraction of trials active for all pyramidal/gad cells

%frac_trials_active = ...
%    cellfun(@(x) cellfun(@(y) nanmean(permute(any(y,2),[1 3 2]),1),x,'uni',false),aligned_df_all,'uni',false);

% only look at some type of trials
% SET TRIAL TYPE HERE
trial_type = 1; % order as before: cl,cr,il,ir
frac_trials_active = ...
    cellfun(@(df,bhv) cellfun(@(curr_df,curr_bhv) ...
    nanmean(permute(any(curr_df(curr_bhv(:,trial_type),:,:),2),[1 3 2]),1),df, ...
    bhv,'uni',false),aligned_df_all,bhv_matrix,'uni',false);


pyr_cells = cell(size(animals));
gad_cells = cell(size(animals));
for curr_animal = 1:length(animals)
    curr_gad = false(size(data_all(curr_animal).im(2).roi_trace_df{1},1),1);
    curr_gad(mice(mice_idx{curr_animal}).interneurons) = true;
    curr_pyr = ~curr_gad;
    pyr_cells{curr_animal} = curr_pyr;
    gad_cells{curr_animal} = curr_gad;
end

frac_trials_active_pyr = cellfun(@(x,y) cellfun(@(z) z(y),x,'uni',false), ...
    frac_trials_active,pyr_cells,'uni',false);
frac_trials_active_gad = cellfun(@(x,y) cellfun(@(z) z(y),x,'uni',false), ...
    frac_trials_active,gad_cells,'uni',false);

avg_frac_trials_active_pyr = cellfun(@(x) nanmean(vertcat(x{:}),2),frac_trials_active_pyr,'uni',false);
max_sessions = max(cellfun(@length,avg_frac_trials_active_pyr));
avg_frac_trials_active_pyr_pad = cellfun(@(x) ...
    padarray(x,[max_sessions-length(x) 0],NaN,'post'),avg_frac_trials_active_pyr,'uni',false);
avg_frac_trials_active_pyr_cat = horzcat(avg_frac_trials_active_pyr_pad{:});

avg_frac_trials_active_gad = cellfun(@(x) nanmean(vertcat(x{:}),2),frac_trials_active_gad,'uni',false);
max_sessions = max(cellfun(@length,avg_frac_trials_active_gad));
avg_frac_trials_active_gad_pad = cellfun(@(x) ...
    padarray(x,[max_sessions-length(x) 0],NaN,'post'),avg_frac_trials_active_gad,'uni',false);
avg_frac_trials_active_gad_cat = horzcat(avg_frac_trials_active_gad_pad{:});



trial_type = 1; % order as before: cl,cr,il,ir
binary_trials_active = ...
    cellfun(@(df,bhv) cellfun(@(curr_df,curr_bhv) ...
    permute(any(curr_df(curr_bhv(:,trial_type),:,:),2),[1 3 2]),df, ...
    bhv,'uni',false),aligned_df_all,bhv_matrix,'uni',false);


%% Activity across reversals

% split aligned activity into pyr/gad
aligned_df_gad = cellfun(@(x,y) cellfun(@(z) z(:,:,mice(y).interneurons),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);
aligned_df_pyr = cellfun(@(x,y) cellfun(@(z) z(:,:,setdiff(1:size(z,3),mice(y).interneurons)),...
    x(cellfun(@(r) ~isempty(r),x)),'uni',false), aligned_df_all,mice_idx,'uni',false);

% get CL activity surrounding reversals (skip ones that are early
% reversals, these were done by mistake)
surrounding_trials = 20;
all_reversals = cellfun(@(x) cellfun(@(y) find(diff(y(:,2)) ~= 0),x,'uni',false), ...
    odor_trials_all,'uni',false); 
reversals = cellfun(@(x) cellfun(@(y) num2cell(y(y > surrounding_trials)), ...
    x,'uni',false),all_reversals,'uni',false); 

% %%% TEMPORARY: pick random reversals
% all_reversals = cellfun(@(x) cellfun(@(y) randi([1 size(y,1)],sum(diff(y(:,2)) ~= 0),1),x,'uni',false), ...
%     odor_trials_all,'uni',false); 
% reversals = cellfun(@(x) cellfun(@(y) num2cell(y(y > surrounding_trials)), ...
%     x,'uni',false),all_reversals,'uni',false); 
% %%%

% %%% TEMPORARY: randomize df
% aligned_df_pyr_rand = cellfun(@(df) cellfun(@(df) ...
%     df(randperm(size(df,1)),:,:),df,'uni',false),aligned_df_pyr,'uni',false);

pre_reversal_activity = cellfun(@(df,bhv,reversal) ...
    cellfun(@(df,bhv,reversal) ...
    cellfun(@(reversal) df(find(bhv(1:reversal-1,1),surrounding_trials, ...
    'last'),:,:),reversal,'uni',false),df,bhv,reversal,'uni',false), ...
    aligned_df_pyr,condition_trials_all,reversals,'uni',false);

post_reversal_activity = cellfun(@(df,bhv,reversal) ...
    cellfun(@(df,bhv,reversal) ...
    cellfun(@(reversal) df(find(bhv(reversal:end,1),surrounding_trials, ...
    'first')+reversal-1,:,:),reversal,'uni',false),df,bhv,reversal,'uni',false), ...
    aligned_df_pyr,condition_trials_all,reversals,'uni',false);

% pull out fraction active trials per cell (only if there are sufficient
% trials pre/post reversal)
use_reversals = cellfun(@(predf,postdf) cellfun(@(predf,postdf) ...
    cellfun(@(predf,postdf) size(predf,1) == surrounding_trials & ...
    size(postdf,1) == surrounding_trials,predf,postdf),predf,postdf,'uni',false), ...
    pre_reversal_activity,post_reversal_activity,'uni',false);

pre_reversal_fractrials = cellfun(@(df,rev) cellfun(@(df,rev) cellfun(@(df)...
    nanmean(permute(any(df,2),[1 3 2])),df(rev), ...
    'uni',false),df,rev,'uni',false),pre_reversal_activity,use_reversals,'uni',false);

post_reversal_fractrials = cellfun(@(df,rev) cellfun(@(df,rev) cellfun(@(df)...
    nanmean(permute(any(df,2),[1 3 2])),df(rev), ...
    'uni',false),df,rev,'uni',false),post_reversal_activity,use_reversals,'uni',false);


% plot frac trials active around reversals
bin_edges = [0:0.2:0.8 Inf];
bin_plot = [0.1:0.2:0.9];
figure; hold on;
line([0 1],[0 1]);
for curr_animal = [1 2 3 6];
    % bin data for errorbar'ing
    if ~isempty(pre_reversal_fractrials{curr_animal}{6})
        [n bins] = histc(pre_reversal_fractrials{curr_animal}{6}{1},bin_edges);
        post_grp = grpstats(post_reversal_fractrials{curr_animal}{6}{1},bins);
        plot(bin_plot(unique(bins)),post_grp,'k');
        
    end
    if ~isempty(pre_reversal_fractrials{curr_animal}{10})
        [n bins] = histc(pre_reversal_fractrials{curr_animal}{10}{1},bin_edges);
        post_grp = grpstats(post_reversal_fractrials{curr_animal}{10}{1},bins);
        plot(bin_plot(unique(bins)),post_grp,'r');
    end
    
    curr_sessions = find(cellfun(@(x) ~isempty(x),pre_reversal_fractrials{curr_animal}));
    for curr_session = curr_sessions(curr_sessions > 10)
        for curr_reversal = length(pre_reversal_fractrials{curr_animal}{curr_session})
            [n bins] = histc(pre_reversal_fractrials{curr_animal} ...
                {curr_session}{curr_reversal},bin_edges);
            post_grp = grpstats(post_reversal_fractrials{curr_animal} ...
                {curr_session}{curr_reversal},bins);
            plot(bin_plot(unique(bins)),post_grp,'g');
        end
    end
end

% plot frac trials active around reversals (mean / errorbar'd)
bin_edges = [0:0.2:0.8 Inf];
bin_plot = [0.1:0.2:0.9];
use_animals = [1 2 3 6];
post_reversal_1 = nan(length(use_animals),length(bin_edges)-1);
post_reversal_2 = nan(length(use_animals),length(bin_edges)-1);
post_reversal_n = nan(length(use_animals),length(bin_edges)-1);
for curr_animal_idx = 1:length(use_animals);
    curr_animal = use_animals(curr_animal_idx);
    % bin data for errorbar'ing
    if ~isempty(pre_reversal_fractrials{curr_animal}{6})
        [n bins] = histc(pre_reversal_fractrials{curr_animal}{6}{1},bin_edges);
        post_grp = grpstats(post_reversal_fractrials{curr_animal}{6}{1},bins);
        post_reversal_1(curr_animal_idx,unique(bins)) = post_grp;        
    end
    
    if ~isempty(pre_reversal_fractrials{curr_animal}{10})
        [n bins] = histc(pre_reversal_fractrials{curr_animal}{10}{1},bin_edges);
        post_grp = grpstats(post_reversal_fractrials{curr_animal}{10}{1},bins);
        post_reversal_2(curr_animal_idx,unique(bins)) = post_grp;
    end
    
    curr_sessions = find(cellfun(@(x) ~isempty(x),pre_reversal_fractrials{curr_animal}));
    curr_reversal_idx = 1;
    curr_post_reversal_n = nan(1,length(bin_edges)-1);
    for curr_session = curr_sessions(curr_sessions > 10)
        for curr_reversal = length(pre_reversal_fractrials{curr_animal}{curr_session})
            [n bins] = histc(pre_reversal_fractrials{curr_animal} ...
                {curr_session}{curr_reversal},bin_edges);
            post_grp = grpstats(post_reversal_fractrials{curr_animal} ...
                {curr_session}{curr_reversal},bins);
            curr_post_reversal_n(curr_reversal_idx,unique(bins)) = post_grp;
            curr_reversal_idx = curr_reversal_idx + 1;
        end
    end
    post_reversal_n(curr_animal_idx,unique(bins)) = nanmean(curr_post_reversal_n);
end

figure; hold on;
line([0 1],[0 1]);
errorbar(bin_plot,nanmean(post_reversal_1),nanstd(post_reversal_1)./sqrt(sum(~isnan(post_reversal_1))),'k')
errorbar(bin_plot,nanmean(post_reversal_2),nanstd(post_reversal_2)./sqrt(sum(~isnan(post_reversal_2))),'r')
errorbar(bin_plot,nanmean(post_reversal_n),nanstd(post_reversal_n)./sqrt(sum(~isnan(post_reversal_n))),'g')

% get general activity across trial
trial_activity = cellfun(@(df,bhv) cellfun(@(df,bhv) permute(any(df(bhv(:,1),:,:),2),[3 1 2]), ...
    df,bhv,'uni',false),aligned_df_pyr,condition_trials_all,'uni',false);





%% Get population changes across reversals

% get general activity across trial
trial_activity = cellfun(@(df,bhv) cellfun(@(df,bhv) permute(any(df(bhv(:,1),:,:),2),[3 1 2]), ...
    df,bhv,'uni',false),aligned_df_pyr,condition_trials_all,'uni',false);

% plot trial-trial active population correlation across reversals
figure;
rev_animals = [1 2 3 6];
curr_sessions = [6 10 11 12];

for cs = 1:length(curr_sessions)
    
    curr_session = curr_sessions(cs);
    col = lines(length(rev_animals));
    smooth_trials = 10;
    smooth_conv = ones(1,smooth_trials)./smooth_trials;
    
    % plot pop corr
    subplot(4,2,(cs*2)-1); hold on;
    for i = 1:length(rev_animals)
        if length(trial_activity{rev_animals(i)}) < curr_session
            continue
        end
        
        curr_rev = find(diff(odor_trials_all{rev_animals(i)}{curr_session}( ...
            condition_trials_all{rev_animals(i)}{curr_session}(:,1),2)) ~= 0,1);
        
        curr_corr = corrcoef(trial_activity{rev_animals(i)}{curr_session});
        curr_corr_mean_prerev = conv(nanmean(curr_corr(:,1:curr_rev),1), ...
            smooth_conv,'valid');
        curr_corr_mean_postrev = conv(nanmean(curr_corr(:,curr_rev+1:end),1), ...
            smooth_conv,'valid');
        
        curr_popcorr_plot = [curr_corr_mean_prerev NaN curr_corr_mean_postrev];
        
        plot((1:length(curr_popcorr_plot)) - curr_rev+smooth_trials-1, ...
            curr_popcorr_plot,'color',col(i,:));
    end
    line([0 0],ylim,'color','k','linestyle','--');
    
    % plot licking corr
    subplot(4,2,(cs*2)); hold on;
    for i = 1:length(rev_animals)
        if length(trial_activity{rev_animals(i)}) < curr_session
            continue
        end
        
        curr_rev = find(diff(odor_trials_all{rev_animals(i)}{curr_session}( ...
            condition_trials_all{rev_animals(i)}{curr_session}(:,1),2)) ~= 0,1);
        
        curr_corr = corrcoef(lick_rate{rev_animals(i)}{curr_session}( ...
            condition_trials_all{rev_animals(i)}{curr_session}(:,1),:)');
        curr_corr_mean_prerev = conv(nanmean(curr_corr(:,1:curr_rev),1), ...
            smooth_conv,'valid');
        curr_corr_mean_postrev = conv(nanmean(curr_corr(:,curr_rev+1:end),1), ...
            smooth_conv,'valid');
        
        curr_lickcorr_plot = [curr_corr_mean_prerev NaN curr_corr_mean_postrev];
        
        plot((1:length(curr_lickcorr_plot)) - curr_rev+smooth_trials-1, ...
            curr_lickcorr_plot,'color',col(i,:));
    end
    line([0 0],ylim,'color','k','linestyle','--');
end

























