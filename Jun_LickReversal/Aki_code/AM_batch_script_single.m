function AM_batch_script_single(m)
    mousename = m.name;
    disp(['Analysing ' mousename '.']);
    updated = AM_batch_get_dFoF_from_img(mousename);
    if(updated || ~exist([mousename '_aligned.mat'],'file'))
        disp(['Analysis files for ' mousename ' is updated.']);
        sessions = AM_batch_behaviordispatcher(mousename);

        %% mouse specific code
        if(strcmp(mousename, 'JL043'))
            N_session = length(sessions);
            for n_session = 1:N_session
                if(sessions(n_session).datenum == 735158)
                    disp('Eliminating the first image from day 2 of JL043.')
                    sessions(n_session).imaged_trials(sessions(n_session).xsg_trials == 1)=0;
                    sessions(n_session).xsg_trials(sessions(n_session).xsg_trials>0) ...
                        = sessions(n_session).xsg_trials(sessions(n_session).xsg_trials>0)-1;
                end
            end
        end
        %%

        sessions = AM_batch_load_behavior(sessions);
        save([mousename '_behavior.mat'], 'sessions');

        sessions = AM_batch_load_analysis(sessions);
        save([mousename '_data.mat'], 'sessions');

        data = AM_align_to_odor_onset(sessions);
        save([mousename '_aligned.mat'], '-struct', 'data');
    else
        disp(['Analysis files for ' mousename ' is NOT updated. Skipping ' mousename '.']);
    end

    aligned_data_file = [mousename '_aligned.mat'];
    vars = whos('-file',aligned_data_file);

    if(~ismember('result',{vars.name}) && ismember('data',{vars.name}))
        disp(['Loading ' aligned_data_file]);
        load(aligned_data_file,'data');
        disp(['Resaving ' aligned_data_file]);
        save(aligned_data_file,'-struct','data');
    end


    if(~ismember('N_trials', {vars.name}) && ismember('N_trial',{vars.name}))
        disp(['Correcting N_trial ' aligned_data_file]);
        load(aligned_data_file,'N_trial');
        N_trials = N_trial;
        N_trial = sum(N_trial);
        save(aligned_data_file,'N_trial','N_trials','-append');
    end

    if(~ismember('hmmstates', {vars.name}) && ismember('result',{vars.name}))
        disp(['Adding HMM states to ' aligned_data_file]);
        load(aligned_data_file,'result');
        hmmstates =  AM_get_hmm_states(result);
        save(aligned_data_file,'hmmstates','-append');
    end

    if(~ismember('first_lick_frame', {vars.name}) && ismember('bhv_frames',{vars.name}))
        disp(['Adding licks to ' aligned_data_file]);
        lick_interval_threshold = 1.2; % frames
        load(aligned_data_file,'bhv_frames');
        load(aligned_data_file,'correct_lick');
        load(aligned_data_file,'imaged_trials');
        data =  AM_get_first_last_lick_frame(bhv_frames,correct_lick>0 & imaged_trials>0, lick_interval_threshold);
        save(aligned_data_file,'-struct','data','-append');
    end
end