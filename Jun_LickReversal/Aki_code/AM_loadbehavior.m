function experiments = AM_loadbehavior(experiments)
    iti_threshold = 10;

    L = size(experiments);
    for i = 1:L
       disp(['Loading ' experiments(i).behavior_file]);
        warning off;
        bhv_data = load(experiments(i).behavior_file);
        warning on;
        experiemnts(i).bhv_data = bhv_data;

        W = length(experiments(i).bhv_frames);  % number of trials
        experiments(i).W = W;

        for w = 1:W
            experiments(i).iti(w) = bhv_data.saved_history.TimesSection_iti{w};
            experiments(i).iti_threshold = iti_threshold;
            experiments(i).iti_ok(w) = bhv_data.saved_history.TimesSection_iti{w} >= iti_threshold;

            experiments(i).applied_odor(w) = bhv_data.saved_history.OdorChoiceSection_AppliedOdor{w};

            experiments(i).rewarded_odor(w) = bhv_data.saved_history.OdorChoiceSection_RewardedOdor{w};

            experiments(i).correct_lick(w) = isfield(experiments(i).bhv_frames{w}.states,'correct_lick') ...
                &&  ~isempty(experiments(i).bhv_frames{w}.states.correct_lick);
            experiments(i).incorrect_rejection(w)  ...
                =experiments(i).applied_odor(w) ==  experiments(i).rewarded_odor(w) ...
                & ~experiments(i).correct_lick(w);


            experiments(i).incorrect_lick(w) = isfield(experiments(i).bhv_frames{w}.states,'incorrect_lick') ...
                &&  ~isempty(experiments(i).bhv_frames{w}.states.incorrect_lick);
            experiments(i).correct_rejection(w)  ...
                =experiments(i).applied_odor(w) ~=  experiments(i).rewarded_odor(w) ...
                & ~experiments(i).incorrect_lick(w);


        end 
    end
    if(nargout == 0)
        save([experiments(1).mousename '_bhv'], 'experiments');
    end
end