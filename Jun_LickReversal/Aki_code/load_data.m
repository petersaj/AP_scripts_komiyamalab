function m = load_data()
    m = mousedata();
    is_valid = false(size(m));
    
    for n_mouse = 1:length(m)
        try
            disp(['Loading ' m(n_mouse).name]);
            data = load([m(n_mouse).name '_aligned.mat']);
            m(n_mouse).data = data;
            if(length(m(n_mouse).data.N_trial) == 1)
                m(n_mouse).data.N_trials = m(n_mouse).data.N_trial;
                m(n_mouse).data.N_trial = sum(m(n_mouse).data.N_trials);
            end
            m(n_mouse).data.N_roi = size(m(n_mouse).data.traces,3);
            m(n_mouse).data.imaged_trials = m(n_mouse).data.imaged_trials & ~cellfun(@isempty,m(n_mouse).data.bhv_frames);
            is_valid(n_mouse) = true;
        catch e
            disp(e.message)
        end
    end
    
%     m = m(is_valid);
end





% 
% 
%     %%
%     for n_mouse = 1:length(m)
%         if(~isempty(m(n_mouse).data))
%             if(length(m(n_mouse).data.N_trial) == 1)
%                 m(n_mouse).data.N_trials = m(n_mouse).data.N_trial;
%                 m(n_mouse).data.N_trial = sum(m(n_mouse).data.N_trials);
%             end
%             m(n_mouse).data.N_roi = size(m(n_mouse).data.traces,3);
%             
% %             m(n_mouse).data.lick_begin = NaN(1,sum(m(n_mouse).data.N_trial));
% %             m(n_mouse).data.lick_end = NaN(1,sum(m(n_mouse).data.N_trial));
% %             
%             N_trial = m(n_mouse).data.N_trial;
% 
%             m(n_mouse).data.lick_frames = cell(N_trial,1);
%             m(n_mouse).data.licks = cell(N_trial,1);
%             m(n_mouse).data.lick_durations = cell(N_trial,1);
%             m(n_mouse).data.intervals_before_lick = cell(N_trial,1);
%             m(n_mouse).data.N_lick = NaN(1,N_trial);
%             m(n_mouse).data.reward_lick = NaN(1,N_trial);
% 
%             m(n_mouse).data.first_lick_frame = NaN(1,N_trial);
%             m(n_mouse).data.reward_lick_frame = NaN(1,N_trial);
%             m(n_mouse).data.last_lick_frame = NaN(1,N_trial);
%             for n_trial = 1:sum(m(n_mouse).data.N_trial) 
%                 bhv_frame = m(n_mouse).data.bhv_frames{n_trial};
%                 if(m(n_mouse).data.correct_lick(n_trial) && m(n_mouse).data.imaged_trials(n_trial))
%                     odor_frame = bhv_frame.states.apply_odor(1);
%                     reward_lick_frame = bhv_frame.states.correct_lick(1);
% 
%                     licks = bhv_frame.pokes.C;
%                     N_lick = size(licks,1);
%                     reward_lick = find(abs(licks-reward_lick_frame)<0.001,1);
% 
% 
%                     lick_duration = diff(licks,1,2);
%                     interval_before_lick = [NaN;diff(licks(:,2))] - lick_duration;
% 
%                     m(n_mouse).data.lick_frames{n_trial} = floor(licks(:,1)-odor_frame)+20;
%                     m(n_mouse).data.N_lick(n_trial) = N_lick;
%                     m(n_mouse).data.reward_lick(n_trial) = reward_lick;
%                     m(n_mouse).data.licks{n_trial} = licks;
%                     m(n_mouse).data.lick_durations{n_trial} = lick_duration;
%                     m(n_mouse).data.intervals_before_lick{n_trial} = interval_before_lick;
%                     
%                     
%                     i = 0;
%                     while(reward_lick > i && m(n_mouse).data.intervals_before_lick{n_trial}(reward_lick-i)<lick_interval_threshold)
%                         i = i+1;
%                     end
%                     first_lick = reward_lick - i;
% 
%                     i = 1;
%                     while(reward_lick + i < N_lick && m(n_mouse).data.intervals_before_lick{n_trial}(reward_lick+i)<lick_interval_threshold)
%                         i = i+1;
%                     end
%                     i = i-1;
%                     last_lick = reward_lick + i;
% 
%                     m(n_mouse).data.first_lick_frame(n_trial) = m(n_mouse).data.lick_frames{n_trial}(first_lick);
%                     m(n_mouse).data.reward_lick_frame(n_trial) = m(n_mouse).data.lick_frames{n_trial}(reward_lick);
%                     m(n_mouse).data.last_lick_frame(n_trial) = m(n_mouse).data.lick_frames{n_trial}(last_lick);
% 
%                 end
%                 
% %                 if(isnan(m(n_mouse).data.lick_frame(1,n_trial)))
% %                     continue;
% %                 end
% % 
% %                 tmp =  find(1==m(n_mouse).data.lick_frame(m(n_mouse).data.odor_frame:end,n_trial),1);
% %                 if(isempty(tmp))
% %                     m(n_mouse).data.lick_begin(n_trial) = NaN;
% %                 else
% %                     m(n_mouse).data.lick_begin(n_trial) = m(n_mouse).data.odor_frame + tmp;
% % 
% %                     tmp = find( 0==m(n_mouse).data.lick_frame((m(n_mouse).data.lick_begin(n_trial)):(end-1),n_trial) ...
% %                               & 0==m(n_mouse).data.lick_frame((m(n_mouse).data.lick_begin(n_trial)+1):end,n_trial) ...
% %                               ,1);
% %                     if(isempty(tmp))
% %                         m(n_mouse).data.lick_end(n_trial) = NaN;
% %                     else
% %                         m(n_mouse).data.lick_end(n_trial) = m(n_mouse).data.lick_begin(n_trial) + tmp;
% %                     end
% %                 end
% 
%             end
%             
%         end
%     end