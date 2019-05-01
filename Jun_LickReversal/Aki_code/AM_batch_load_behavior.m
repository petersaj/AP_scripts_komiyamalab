function sessions = AM_batch_load_behavior(sessions)
    iti_threshold = 10;

    L = size(sessions);
    for i = 1:L
       disp(['Loading ' sessions(i).behavior_file]);
        warning off;
        bhv_data = load(sessions(i).behavior_file);
        warning on;
        sessions(i).bhv_data = bhv_data;

        W = length(sessions(i).bhv_frames);  % number of trials
        sessions(i).W = W;
        
        sessions(i).lick_frame = [];
        sessions(i).lick_length = [];
        sessions(i).odor_frame = [];
        sessions(i).reward_frame = [];
        sessions(i).lick_after_reward_frame = [];
        
        for w = 1:W
            sessions(i).iti(w) = bhv_data.saved_history.TimesSection_iti{w};
            sessions(i).iti_threshold = iti_threshold;
            sessions(i).iti_ok(w) = bhv_data.saved_history.TimesSection_iti{w} >= iti_threshold;

            sessions(i).applied_odor(w) = bhv_data.saved_history.OdorChoiceSection_AppliedOdor{w};

            sessions(i).rewarded_odor(w) = bhv_data.saved_history.OdorChoiceSection_RewardedOdor{w};

            sessions(i).correct_lick(w) = isfield(sessions(i).bhv_frames{w}.states,'correct_lick') ...
                &&  ~isempty(sessions(i).bhv_frames{w}.states.correct_lick);
            sessions(i).incorrect_rejection(w)  ...
                =sessions(i).applied_odor(w) ==  sessions(i).rewarded_odor(w) ...
                & ~sessions(i).correct_lick(w);
            
            
            sessions(i).incorrect_lick(w) = isfield(sessions(i).bhv_frames{w}.states,'incorrect_lick') ...
                &&  ~isempty(sessions(i).bhv_frames{w}.states.incorrect_lick);
            sessions(i).correct_rejection(w)  ...
                =sessions(i).applied_odor(w) ~=  sessions(i).rewarded_odor(w) ...
                & ~sessions(i).incorrect_lick(w);
            
            if(sessions(i).imaged_trials(w))
                sessions(i).lick_frame = [sessions(i).lick_frame; floor(sessions(i).bhv_frames{w}.pokes.C...
                    (~isnan(sessions(i).bhv_frames{w}.pokes.C(:,1)) ,1))+(sessions(i).xsg_trials(w)-1)*1000];
                
                sessions(i).lick_length = [sessions(i).lick_length; ...
                    sessions(i).bhv_frames{w}.pokes.C(~isnan(sessions(i).bhv_frames{w}.pokes.C(:,1)),2)-...
                    sessions(i).bhv_frames{w}.pokes.C(~isnan(sessions(i).bhv_frames{w}.pokes.C(:,1)),1)];
                
                sessions(i).odor_frame = [sessions(i).odor_frame; floor(sessions(i).bhv_frames{w}.states.apply_odor(1))+...
                    (sessions(i).xsg_trials(w)-1)*1000];
                if(isfield(sessions(i).bhv_frames{w}.states,'correct_lick')&&~isempty(sessions(i).bhv_frames{w}.states.correct_lick))
                    sessions(i).reward_frame = [sessions(i).reward_frame; floor(sessions(i).bhv_frames{w}.states.correct_lick(1))+...
                        (sessions(i).xsg_trials(w)-1)*1000];
                    
                    sessions(i).lick_after_reward_frame = [sessions(i).lick_after_reward_frame; floor(sessions(i).bhv_frames{w}.pokes.C...
                        (find(sessions(i).bhv_frames{w}.pokes.C(:,1) > sessions(i).bhv_frames{w}.states.correct_lick(1)+0.001,1),1))+(sessions(i).xsg_trials(w)-1)*1000];

                end
                
                sessions(i).trial_num( (floor(sessions(i).bhv_frames{w}.states.state_0(3))...
                    :floor(sessions(i).bhv_frames{w}.states.state_0(2)))+(sessions(i).xsg_trials(w)-1)*1000  ) = w;
            else
                  sessions(i).odor_frame = [sessions(i).odor_frame; NaN];
                
            end
            
        end 
    end
    if(nargout == 0)
        save([sessions(1).mousename '_bhv'], 'sessions');
    end
end