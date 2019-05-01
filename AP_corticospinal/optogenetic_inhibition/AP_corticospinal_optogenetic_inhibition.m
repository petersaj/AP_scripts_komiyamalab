%% Load behavior data from optogenetic inhibition animals

animals = {'AP164','AP165','AP166','AP167','AP168', ...
    'AP170','AP171','AP172','AP173'};

% Define whether opto days in each animal were blocked (kwik-cast)
block_days = { ...
    [0,1]; ... % AP164
    [1,0]; ... % AP165
    [0,1]; ... % AP166
    [1,0]; ... % AP167
    [0,1]; ... % AP168
    [0,1]; ... % AP170
    [1,0]; ... % AP171
    [0,1]; ... % AP172
    [1,0]}';    % AP173

block_days = cellfun(@(x) logical(x),block_days,'uni',false);

analysis = AP_corticospinal_process_bhv_only(animals);

% Only use days with opto trials
for curr_animal = 1:length(analysis)
    nonopto_sessions = cellfun(@(x) isempty(x),{analysis(curr_animal).lever.opto_trial});
    analysis(curr_animal).lever(nonopto_sessions) = [];
    analysis(curr_animal).bhv(nonopto_sessions) = [];
end


figure; 

% Fraction of rewarded trials for opto vs. non-opto

clear_rewarded_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).rewarded_trials),{analysis.bhv},block_days,'uni',false);
clear_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).opto_trials),{analysis.bhv},block_days,'uni',false);

blocked_rewarded_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).rewarded_trials),{analysis.bhv},block_days,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).opto_trials),{analysis.bhv},block_days,'uni',false);

clear_rewarded_opto_mean = cellfun(@(rew,opto) ...
    nanmean(rew(opto),1),clear_rewarded_trials_cat,clear_opto_trials_cat);
clear_rewarded_nonopto_mean = cellfun(@(rew,opto) ...
    nanmean(rew(~opto),1),clear_rewarded_trials_cat,clear_opto_trials_cat);

blocked_rewarded_opto_mean = cellfun(@(rew,opto) ...
    nanmean(rew(opto),1),blocked_rewarded_trials_cat,blocked_opto_trials_cat);
blocked_rewarded_nonopto_mean = cellfun(@(rew,opto) ...
    nanmean(rew(~opto),1),blocked_rewarded_trials_cat,blocked_opto_trials_cat);

subplot(5,3,1); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_nonopto_mean;clear_rewarded_opto_mean]);
plot([clear_rewarded_nonopto_mean;clear_rewarded_opto_mean],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Fraction correct trials')
title('Clear window')

subplot(5,3,2); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([blocked_rewarded_nonopto_mean;blocked_rewarded_opto_mean]);
plot([blocked_rewarded_nonopto_mean;blocked_rewarded_opto_mean],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Fraction correct trials')
title('Blocked window')

subplot(5,3,3); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_opto_mean-clear_rewarded_nonopto_mean; ...
    blocked_rewarded_opto_mean-blocked_rewarded_nonopto_mean]);
plot([clear_rewarded_opto_mean-clear_rewarded_nonopto_mean; ...
    blocked_rewarded_opto_mean-blocked_rewarded_nonopto_mean],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Clear','Blocked'});
ylabel('\Delta fraction correct trials')
title('Effect of light')


% Cue to reward time

clear_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).cue_to_reward_time),{analysis.lever},block_days,'uni',false);
clear_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).opto_trial),{analysis.lever},block_days,'uni',false);

blocked_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).cue_to_reward_time),{analysis.lever},block_days,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).opto_trial),{analysis.lever},block_days,'uni',false);

clear_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);
blocked_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);

clear_time_cat = cellfun(@(x,cued) x(cued),clear_time_cat,clear_cued_trials_cat,'uni',false);
clear_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    clear_opto_trials_cat,clear_cued_trials_cat,'uni',false);
blocked_time_cat = cellfun(@(x,cued) x(cued),blocked_time_cat,blocked_cued_trials_cat,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    blocked_opto_trials_cat,blocked_cued_trials_cat,'uni',false);

clear_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),clear_time_cat,clear_opto_trials_cat);
clear_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),clear_time_cat,clear_opto_trials_cat);

blocked_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),blocked_time_cat,blocked_opto_trials_cat);
blocked_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),blocked_time_cat,blocked_opto_trials_cat);

subplot(5,3,4); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median]);
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Cue to reward time')
title('Clear window')

subplot(5,3,5); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median]);
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Cue to reward time')
title('Blocked window')

subplot(5,3,6); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median]);
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Clear','Blocked'});
ylabel('\Delta cue to reward time')
title('Effect of light')

% Reaction time

clear_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).reaction_time),{analysis.lever},block_days,'uni',false);
clear_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).opto_trial),{analysis.lever},block_days,'uni',false);

blocked_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).reaction_time),{analysis.lever},block_days,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).opto_trial),{analysis.lever},block_days,'uni',false);

clear_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);
blocked_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);

clear_time_cat = cellfun(@(x,cued) x(cued),clear_time_cat,clear_cued_trials_cat,'uni',false);
clear_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    clear_opto_trials_cat,clear_cued_trials_cat,'uni',false);
blocked_time_cat = cellfun(@(x,cued) x(cued),blocked_time_cat,blocked_cued_trials_cat,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    blocked_opto_trials_cat,blocked_cued_trials_cat,'uni',false);

clear_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),clear_time_cat,clear_opto_trials_cat);
clear_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),clear_time_cat,clear_opto_trials_cat);

blocked_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),blocked_time_cat,blocked_opto_trials_cat);
blocked_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),blocked_time_cat,blocked_opto_trials_cat);

subplot(5,3,7); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median]);
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Reaction time')
title('Clear window')

subplot(5,3,8); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median]);
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Reaction time')
title('Blocked window')

subplot(5,3,9); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median]);
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Clear','Blocked'});
ylabel('\Delta reaction time')
title('Effect of light')


% Move onset to reward time

clear_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).move_to_reward_time),{analysis.lever},block_days,'uni',false);
clear_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).opto_trial),{analysis.lever},block_days,'uni',false);

blocked_time_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).move_to_reward_time),{analysis.lever},block_days,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).opto_trial),{analysis.lever},block_days,'uni',false);

clear_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);
blocked_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);

clear_time_cat = cellfun(@(x,cued) x(cued),clear_time_cat,clear_cued_trials_cat,'uni',false);
clear_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    clear_opto_trials_cat,clear_cued_trials_cat,'uni',false);
blocked_time_cat = cellfun(@(x,cued) x(cued),blocked_time_cat,blocked_cued_trials_cat,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    blocked_opto_trials_cat,blocked_cued_trials_cat,'uni',false);

clear_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),clear_time_cat,clear_opto_trials_cat);
clear_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),clear_time_cat,clear_opto_trials_cat);

blocked_rewarded_opto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(opto),1),blocked_time_cat,blocked_opto_trials_cat);
blocked_rewarded_nonopto_median = cellfun(@(rew,opto) ...
    nanmedian(rew(~opto),1),blocked_time_cat,blocked_opto_trials_cat);

subplot(5,3,10); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median]);
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Move to reward time')
title('Clear window')

subplot(5,3,11); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median]);
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('Move to reward time')
title('Blocked window')

subplot(5,3,12); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median]);
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Clear','Blocked'});
ylabel('\Delta move to reward time')
title('Effect of light')


% Movement correlation

clear_move_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).rewarded_movement),{analysis.lever},block_days,'uni',false);
clear_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).opto_trial),{analysis.lever},block_days,'uni',false);

blocked_move_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).rewarded_movement),{analysis.lever},block_days,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).opto_trial),{analysis.lever},block_days,'uni',false);

clear_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(~blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);
blocked_cued_trials_cat = cellfun(@(x,blocked) ...
    vertcat(x(blocked).cued_movement_trials),{analysis.lever},block_days,'uni',false);

clear_move_cat = cellfun(@(x,cued) x(cued),clear_move_cat,clear_cued_trials_cat,'uni',false);
clear_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    clear_opto_trials_cat,clear_cued_trials_cat,'uni',false);
blocked_move_cat = cellfun(@(x,cued) x(cued),blocked_move_cat,blocked_cued_trials_cat,'uni',false);
blocked_opto_trials_cat = cellfun(@(x,cued) logical(x(cued)), ...
    blocked_opto_trials_cat,blocked_cued_trials_cat,'uni',false);

clear_move_cat_cat = cellfun(@(x) horzcat(x{:}),clear_move_cat,'uni',false);
blocked_move_cat_cat = cellfun(@(x) horzcat(x{:}),blocked_move_cat,'uni',false);

clear_rewarded_opto_median = cellfun(@(rew,opto) nanmedian(AP_itril(corrcoef(horzcat( ...
    rew(501:end,opto))),-1)),clear_move_cat_cat,clear_opto_trials_cat);
clear_rewarded_nonopto_median = cellfun(@(rew,opto) nanmedian(AP_itril(corrcoef(horzcat( ...
    rew(501:end,~opto))),-1)),clear_move_cat_cat,clear_opto_trials_cat);

blocked_rewarded_opto_median = cellfun(@(rew,opto) nanmedian(AP_itril(corrcoef(horzcat( ...
    rew(501:end,opto))),-1)),blocked_move_cat_cat,blocked_opto_trials_cat);
blocked_rewarded_nonopto_median = cellfun(@(rew,opto) nanmedian(AP_itril(corrcoef(horzcat( ...
    rew(501:end,~opto))),-1)),blocked_move_cat_cat,blocked_opto_trials_cat);

subplot(5,3,13); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median]);
plot([clear_rewarded_nonopto_median;clear_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('median movement correlation')
title('Clear window')

subplot(5,3,14); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median]);
plot([blocked_rewarded_nonopto_median;blocked_rewarded_opto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'No Light','Light'});
ylabel('median movement correlation')
title('Blocked window')

subplot(5,3,15); hold on;
set(gca,'ColorOrder',jet(length(animals)));
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median]);
plot([clear_rewarded_opto_median-clear_rewarded_nonopto_median; ...
    blocked_rewarded_opto_median-blocked_rewarded_nonopto_median],'o');
xlim([0,3]);
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Clear','Blocked'});
ylabel('\Delta median movement correlation')
title('Effect of light')












