function classified_rois = AP_classify_movement_cells_continuous(data,analysis,split_time)
% classified_rois = AP_classify_movement_cells_continuous(data,analysis,split_time)
%
% Classify movement-related cells using
%
% Inputs: 
% data - from AP_corticospinal_load_all
% analysis - from AP_corticospinal_prepare_loaded
% split_time - can either be done on a session basis (no 3rd input), 
% or within a session split up into time blocks given by split_time
%
% Outputs:
% classified_rois - structure with logical vectors for movement or
% quiescent classifiation and p-values

if nargin < 3
    
    classified_rois = struct('movement',cell(length(data),1),'quiescent',cell(length(data),1), ...
        'p',cell(length(data),1));
    
    for curr_animal = 1:length(data);
        
        n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
        n_sessions = length(data(curr_animal).im);
        
        classified_rois(curr_animal).movement = false(n_rois,n_sessions);
        classified_rois(curr_animal).quiescent = false(n_rois,n_sessions);
        classified_rois(curr_animal).p = nan(n_rois,n_sessions);
        
        for curr_session = 1:n_sessions
                       
            % This frame restriction is only for one old animal
            use_frames = 1:min(length(analysis(curr_animal).lever(curr_session).lever_move_frames), ...
                size(data(curr_animal).im(curr_session).roi_trace_thresh,2));
            
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHORT TIME TEST
%             
%             min_frames = min(cellfun(@(x) size(x,2),{data(curr_animal).im.roi_trace_df}));
%             curr_n_frames = size(data(curr_animal).im(curr_session).roi_trace_df,2);
%             use_frames = 1:curr_n_frames;
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            activity_frames = +(data(curr_animal).im(curr_session).roi_trace_thresh(:,use_frames) > 0);
            
            lever_active_frames = analysis(curr_animal).lever(curr_session).lever_move_frames(use_frames);
            
            % Split active frames by movement blocks / quiescent frames
            % (to split up by movement blocks / quiescent blocks)
            %lever_active_frames_nan = +lever_active_frames;
            %lever_active_frames_nan(~lever_active_frames) = NaN;
            boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
            lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
            
%             %%%%%% TESTING %%%%%%%%%%%%%%%%%%%
%             % restrict movements epochs used to classify by duration
%             
%             % identify short movement epochs
%             move_duration_cutoff = 28*1;
%             frames_split = mat2cell(transpose(1:length(lever_active_frames)),diff(boundary_frames));
%             frame_split_remove = cellfun(@(x) length(x),frames_split) < move_duration_cutoff;
%             split_move = cellfun(@(x) all(x),lever_active_frames_split);
%             
%             cut_split = frame_split_remove & split_move;
%             
%             % set lever times during short movements to 0
%             lever_active_frames(vertcat(frames_split{cut_split})) = 0;
%             lever_active_frames_split(cut_split) = cellfun(@(x) ...
%                 zeros(size(x)), lever_active_frames_split(cut_split),'uni',false);
%             
%             % set activity during short movements to 0
%             activity_frames(:,vertcat(frames_split{cut_split})) = 0;
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
            % Get shuffled activity distribution
            num_rep = 10000;
            shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
            lever_active_shuffle = nan(length(lever_active_frames),num_rep);
            for i = 1:num_rep
                lever_active_shuffle(:,i) = ...
                    vertcat(lever_active_frames_split{shuffle_perms(:,i)});
            end
            clear shuffle_perms
            
            %%% Classify ROIs according to movement/quiescence preference
            movement_activity = activity_frames*lever_active_frames;           
            
            shuffle_movement_activity = activity_frames*lever_active_shuffle;
            clear lever_active_shuffle
            
            movement_rank = tiedrank([movement_activity shuffle_movement_activity]')';
            movement_p = movement_rank(:,1)/(num_rep+1);
            movement_cells = movement_p > 0.975;
            quiescent_cells = movement_p < 0.025;
            
            classified_rois(curr_animal).movement(:,curr_session) = movement_cells;
            classified_rois(curr_animal).quiescent(:,curr_session) = quiescent_cells;            
            classified_rois(curr_animal).p(:,curr_session) = movement_p;
            
            % ROIs that are active but not movement/quiescent classified
            m = classified_rois(curr_animal).movement(:,curr_session);
            q = classified_rois(curr_animal).quiescent(:,curr_session);
            u = ~m & ~q;
            
            m_avg = nanmean(data(curr_animal).im(curr_session).roi_trace_thresh(m,:),2);
            q_avg = nanmean(data(curr_animal).im(curr_session).roi_trace_thresh(q,:),2);
            c_avg = nanmean(data(curr_animal).im(curr_session).roi_trace_thresh(m|q,:),2);
            u_avg = nanmean(data(curr_animal).im(curr_session).roi_trace_thresh(u,:),2);
            
            all_onsets =  sum(diff(data(curr_animal).im(curr_session).roi_trace_thresh > 0,[],2) == 1,2);
            m_onsets = sum(diff(data(curr_animal).im(curr_session).roi_trace_thresh(m,:) > 0,[],2) == 1,2);
            q_onsets = sum(diff(data(curr_animal).im(curr_session).roi_trace_thresh(q,:) > 0,[],2) == 1,2);
            c_onsets = sum(diff(data(curr_animal).im(curr_session).roi_trace_thresh(m|q,:) > 0,[],2) == 1,2);
            u_onsets = sum(diff(data(curr_animal).im(curr_session).roi_trace_thresh(u,:) > 0,[],2) == 1,2);
            
            % Define "active unclassified"
            % 1) using minimum
            %u_active = u_avg > min(c_avg) & u_onsets > min(c_onsets);
            % 2) using percentile
            u_active = u_avg > prctile(c_avg,5) & u_onsets > prctile(c_onsets,5);

            classified_rois(curr_animal).unclassified_active(u,curr_session) = u_active;
                  
            % Secondary criteria: minimum number of discrete events
            min_events = 5;
            classified_rois(curr_animal).movement(all_onsets < min_events,curr_session) = false;
            classified_rois(curr_animal).quiescent(all_onsets < min_events,curr_session) = false;
            classified_rois(curr_animal).unclassified_active(all_onsets < min_events,curr_session) = false;
            
            disp(['Classified session ' num2str(curr_session)]);
        end
        
        disp(['Classified ' data(curr_animal).animal]);
    end
    disp('Finished all')
    
else
    
    classified_rois = struct('movement',cell(length(data),1),'quiescent',cell(length(data),1));
    
    for curr_animal = 1:length(data);
        
        n_rois = size(data(curr_animal).im(1).roi_trace_df,1);
        n_sessions = length(data(curr_animal).im);
        
        classified_rois(curr_animal).movement = cell(n_sessions,1);
        classified_rois(curr_animal).quiescent = cell(n_sessions,1);
        classified_rois(curr_animal).p = cell(n_sessions,1);
        
        for curr_session = 1:n_sessions
            
            [num_rois,num_frames] = size(data(curr_animal).im(curr_session).roi_trace_df);
            
            % Split frames into time bins, overlapping at end to keep time
            frames_split = ceil(linspace(0,num_frames,split_time+1));
            
            act_split = mat2cell(+(data(curr_animal).im(curr_session). ...
                    roi_trace_thresh > 0),num_rois,diff(frames_split));
                
            lever_active_split = mat2cell(+analysis(curr_animal).lever(curr_session). ...
                    lever_move_frames,diff(frames_split),1);
            
            for curr_split = 1:split_time;
                
                curr_act = act_split{curr_split};
                
                lever_active_frames = lever_active_split{curr_split};
                
                % Split active frames by movement blocks / quiescent frames
                % (to split up by movement blocks / quiescent blocks)
                %lever_active_frames_nan = +lever_active_frames;
                %lever_active_frames_nan(~lever_active_frames) = NaN;
                boundary_frames = find(diff([Inf;lever_active_frames;Inf]) ~= 0);
                lever_active_frames_split = mat2cell(lever_active_frames,diff(boundary_frames));
                
                % Get activity during movement/quiescence
                movement_activity = curr_act*lever_active_frames;
                
                % Get shuffled activity distribution
                num_rep = 10000;
                shuffle_perms = shake(repmat([1:length(lever_active_frames_split)]',1,num_rep),1);
                lever_active_shuffle = nan(length(lever_active_frames),num_rep);
                for i = 1:num_rep
                    lever_active_shuffle(:,i) = ...
                        vertcat(lever_active_frames_split{shuffle_perms(:,i)});
                end
                clear shuffle_perms
                
                shuffle_movement_activity = curr_act*lever_active_shuffle;
                clear lever_active_shuffle
                
                movement_rank = tiedrank([movement_activity shuffle_movement_activity]')';
                movement_p = movement_rank(:,1)/(num_rep+1);
                movement_cells = movement_p > 0.975;
                quiescent_cells = movement_p < 0.025;
                
                classified_rois(curr_animal).movement{curr_session}(:,curr_split) = movement_cells;
                classified_rois(curr_animal).quiescent{curr_session}(:,curr_split) = quiescent_cells;                
                classified_rois(curr_animal).p{curr_session}(:,curr_split) = movement_p;
                
            end
            
            disp(['Classified session ' num2str(curr_session)]);
        end
        
        disp(['Classified ' data(curr_animal).animal]);
    end
    disp('Finished all')
    
    
end




















