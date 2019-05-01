% Analyze behavior of prepared data
% data from AP_corticospinal_prepare_processed
% (or AP_corticospinal_process_bhv_only)

num_animals = length(analysis);

% Plot the reaction time
median_reaction_time = arrayfun(@(x) cellfun(@(param,trials) median(param(trials)), ...
    {analysis(x).lever.reaction_time},{analysis(x).lever.cued_movement_trials}),1:num_animals,'uni',false);

% Plot the movement onset to reward time
median_move_to_reward_time = arrayfun(@(x) cellfun(@(param,trials) median(param(trials)), ...
    {analysis(x).lever.move_to_reward_time},{analysis(x).lever.cued_movement_trials}),1:num_animals,'uni',false);

% Plot the movement time
median_movement_time = arrayfun(@(x) cellfun(@(param,trials) median(param(trials)), ...
    {analysis(x).lever.rewarded_movement_time},{analysis(x).lever.cued_movement_trials}),1:num_animals,'uni',false);

% Plot the trial-by-trial movement correlation
movement_start_time = 1001; % (defined in prepare_processed);
movement_use_time = 3000; % ms 

% (do in loop in case dispatcher file not saved, was in old cohort)
max_sessions = max(arrayfun(@(x) length(analysis(x).lever),1:num_animals));
movement_use_trials = cell(num_animals,1);
for curr_animal = 1:num_animals    
   for curr_session = 1:length(analysis(curr_animal).lever)
      if ~isempty(analysis(curr_animal).lever(curr_session).rewarded_movement)
          
          curr_move = ...
              horzcat(analysis(curr_animal).lever(curr_session).rewarded_movement{ ...
              analysis(curr_animal).lever(curr_session).cued_movement_trials});
          
          movement_use_trials{curr_animal}{curr_session} = curr_move;          
      end
   end
end

movement_corr_grid = nan(max_sessions,max_sessions,num_animals);
for curr_animal = 1:num_animals
    for session_x = 1:length(analysis(curr_animal).lever);
        for session_y = 1:length(analysis(curr_animal).lever)            
            if ~isempty(movement_use_trials{curr_animal}{session_x}) && ...
                    ~isempty(movement_use_trials{curr_animal}{session_y})
                
                num_trials_x = size(movement_use_trials{curr_animal}{session_x},2);
                
                curr_move_corr = corrcoef([movement_use_trials{ ...
                    curr_animal}{session_x}(movement_start_time:movement_start_time+movement_use_time,:) ...
                    movement_use_trials{curr_animal}{session_y}( ...
                    movement_start_time:movement_start_time+movement_use_time,:)]);
                
                curr_use_corr = curr_move_corr(1:num_trials_x,num_trials_x+1:end);
                
                movement_corr_grid(session_x,session_y,curr_animal) = nanmedian(curr_use_corr(:));
                
            end
        end
    end
end


%% Plot behavior summary

% Timing (movement time, reaction time, move to reward time)
figure; 

movement_time_cat = vertcat(median_movement_time{:});

subplot(3,2,1)
plot(movement_time_cat','linewidth',2);
xlim([0 15]);
ylabel('Median movement time (s)')
xlabel('Day')
subplot(3,2,2);
errorbar(nanmean(movement_time_cat),nanstd(movement_time_cat)./ ...
    sqrt(sum(~isnan(movement_time_cat))),'k','linewidth',2);
xlim([0 15]);
ylabel('Median movement time (s)')
xlabel('Day')

reaction_time_cat = vertcat(median_reaction_time{:});

subplot(3,2,3)
plot(reaction_time_cat','linewidth',2);
xlim([0 15]);
ylabel('Median reaction time (s)')
xlabel('Day')
subplot(3,2,4);
errorbar(nanmean(reaction_time_cat),nanstd(reaction_time_cat)./ ...
    sqrt(sum(~isnan(reaction_time_cat))),'k','linewidth',2);
xlim([0 15]);
ylabel('Median reaction time (s)')
xlabel('Day')

move_to_reward_time_cat = vertcat(median_move_to_reward_time{:});

subplot(3,2,5)
plot(move_to_reward_time_cat','linewidth',2);
xlim([0 15]);
ylabel('Median move to reward time (s)')
xlabel('Day')
subplot(3,2,6);
errorbar(nanmean(move_to_reward_time_cat),nanstd(move_to_reward_time_cat)./ ...
    sqrt(sum(~isnan(move_to_reward_time_cat))),'k','linewidth',2);
xlim([0 15]);
ylabel('Median move to reward time (s)')
xlabel('Day')


% Movement correlation
figure; 
subplot(1,3,1);
imagesc(nanmean(movement_corr_grid,3));colormap(hot);
ylabel('Day');
xlabel('Day');
title('Movement correlation');

first_movecorr_diag = cell2mat(arrayfun(@(x) diag(movement_corr_grid(:,:,x)), ...
    1:size(movement_corr_grid,3),'uni',false))';
second_movecorr_diag = cell2mat(arrayfun(@(x) diag(movement_corr_grid(:,:,x),-1), ...
    1:size(movement_corr_grid,3),'uni',false))';

subplot(2,3,2);
plot(first_movecorr_diag','linewidth',2);
xlim([0 15]);
ylabel('Median movement correlation')
xlabel('Day');
title('Within day');

subplot(2,3,3);
errorbar(nanmean(first_movecorr_diag),nanstd(first_movecorr_diag)./ ...
    sqrt(sum(~isnan(first_movecorr_diag))),'k','linewidth',2);
xlim([0 15]);
ylabel('Median movement correlation')
xlabel('Day');
title('Within day');

subplot(2,3,5);
plot(second_movecorr_diag','linewidth',2);
xlim([0 14]);
ylabel('Median movement correlation')
xlabel('Day v. next day');
title('Next day');

subplot(2,3,6);
errorbar(nanmean(second_movecorr_diag),nanstd(second_movecorr_diag)./ ...
    sqrt(sum(~isnan(second_movecorr_diag))),'k','linewidth',2);
xlim([0 14]);
ylabel('Median movement correlation')
xlabel('Day v. next day');
title('Next day');
















