function [peak_frames peak_amplitudes] = AP_schmittDetect2(v, delta, peak_back, local_max, drop_percentage)
% function finds peaks that satisfies schmitt trigger criterion, [peak_frames] = ap_schmittDetect(v, delta, peak_back, local_max, drop_percentage)

% v is vector to detect next peak
% delta is minimum difference required, can be vector of any length
% peak_back is how many frames ago the first jump had to occur within
% local_max is finding the peak point within this many points in front
% drop_percentage is percent amplitude it has to drop to before finding new

% Look peak_frames back for a jump in delta(1)
schmitt_matrix = repmat(v,peak_back,1);
for i = 1:size(schmitt_matrix,1)
    schmitt_matrix(i,:) = circshift(schmitt_matrix(i,:),[0 i]);
end
schmitt_matrix = repmat(v,peak_back,1) - schmitt_matrix;
candidate_frames = find([0 diff(any(schmitt_matrix > delta(1))) == 1])+1;

% get rid of uncheckable candidate frames 
candidate_frames(candidate_frames < peak_back+1) = [];
candidate_frames(candidate_frames > length(v)-length(delta)-1-local_max) = [];

peak_frames = [];
peak_amplitudes = [];
curr_amplitude = Inf;
drop_flag = 1;
pre_baseline_frames = 10;

for i = candidate_frames;
    
    % find where vector goes above critera within peak_back frames
    prior_frames = v(i-peak_back:i-1);
    prior_frames_diff = max(v(i:i+local_max)) - prior_frames;
    prior_min = find(prior_frames_diff == max(prior_frames_diff),1);
    prior_min = (peak_back+1) - prior_min; % prior_min gives minimum in frames back

    
%      %%%% FIX LATER: at the moment this doesn't work, false negatives
%     % last minimum can't be >2 stds from avg of all of them: protects from
%     % quick downward artifacts
%     downward_artifact = 0;
%     mean_pre_baseline_frames = mean(v(i-(prior_min+pre_baseline_frames+1):i-prior_min-1));
%     std_pre_baseline_frames = std(v(i-(prior_min+pre_baseline_frames+1):i-prior_min-1));
%     if abs(v(i-prior_min)) > mean_pre_baseline_frames+2*std_pre_baseline_frames
%         downward_artifact = 1;
%     end
    
%     % flag that trace has dipped below drop_percentage of last amplitude,
%     % UNLESS IT DROPS JUST FOR LESS THAN 5 FRAMES: that's prob an artifact,
%     % also if it find the drop in this part then you have to skip over it
%     % to look for the next peak
%     if drop_flag == 0;
%         if all(v(i:i+5) <= drop_percentage*curr_amplitude);
%             drop_flag = 1;
%             i = i+5;
%             continue
%         end
%     end

    if any(prior_frames_diff >= delta(1));
        if all(v(i:i+length(delta)-1) >= delta+v(i-prior_min)) && ...
                i+local_max <= length(v);
            % check on flags
            if drop_flag == 1 %&& downward_artifact == 0;
                
                % the start of the transient is where the jump was < 20%
                peak_start = (peak_back+1)-find(prior_frames <= 0.2,1,'last');
                % if the trace didn't go below 20%, just use minimum
                if isempty(peak_start)
                     peak_start = (peak_back+1)-find(prior_frames ...
                        == min(prior_frames),1,'last');
                end
                peak_frames = [peak_frames i-peak_start];
                curr_amplitude = max(v(i:i+local_max));
                peak_amplitudes = [peak_amplitudes curr_amplitude];
                if length(i-peak_start) ~= length(curr_amplitude)
                    disp('ap_schmittDetect, frames ~= amplitudes!')
                    keyboard
                end
                % set the drop flag back to 1, wait until it drops
                %drop_flag = 0;
                %i = find(v(i:end) == max(v(i:i+local_max)),1)+i+1;
            else
                continue
            end
        else
            continue
        end
    else
        continue
    end    
end

% get rid of duplicated peaks
[~, unique_peak_idx, ~] = unique(peak_frames,'first');
peak_frames = peak_frames(unique_peak_idx);
peak_amplitudes = peak_amplitudes(unique_peak_idx);
