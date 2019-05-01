%% Analysis test script for Jun


aligned_df_all

condition_trials_all



CL_activity = aligned_df_all{mouse}{session}(CL,:,:)
CL_activity_binary = any(CL_activity,2);
CL_activity_binary_mean = mean(CL_activity_binary,1);

















