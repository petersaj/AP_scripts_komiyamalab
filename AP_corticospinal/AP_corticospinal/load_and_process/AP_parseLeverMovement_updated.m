function [lever_active,lever_force_resample,lever_force_smooth, ...
    lever_velocity_envelope_smooth] = AP_parseLeverMovement_updated(xsg_data,downsampled)
%[lever_active,lever_force_resample,lever_force_smooth, ...
%    lever_velocity_envelope_smooth] = AP_parseLeverMovement_updated(xsg_data,downsampled)
%
% parse the lever force into movement/nonmovement epochs
%
% 'xsg_data' can either be the structure from looped data, or raw lever trace
% 'downsampled' is flag if lever is already at 1 kHz

if isstruct(xsg_data)
    % If xsg_data is structure from looped data, pull out sample rate
    xsg_sample_rate = xsg_data.header.acquirer.acquirer.sampleRate;
    lever_force_raw = xsg_data.data.acquirer.trace_2;
else
    % If not (continuous imaging), pull out lever during imaged time
    lever_force_raw = xsg_data;
    xsg_sample_rate = 10000;
end

% resample lever trace to 1kHz
if nargin == 1 || ~downsampled
    [n d] = rat(1000/xsg_sample_rate);
    lever_force_resample = resample(lever_force_raw,n,d);
elseif downsampled
    lever_force_resample = lever_force_raw;
end

butterworth_stop = 5/500; % fraction of nyquist (cutoff = 10 Hz)
[b a] = butter(4, butterworth_stop,'low');
lever_force_smooth = filtfilt(b,a,lever_force_resample);

lever_velocity_resample = [0;diff(lever_force_smooth)];
lever_velocity_resample_smooth = smooth(lever_velocity_resample,5);
lever_velocity_resample_smooth(isnan(lever_velocity_resample)) = NaN;

% get velocity envelope, smooth over 5 frames (~150ms)
lever_velocity_hilbert = hilbert(lever_velocity_resample_smooth);
lever_velocity_envelope = sqrt(lever_velocity_hilbert.* ...
    conj(lever_velocity_hilbert));
% at the moment - don't smooth
lever_velocity_envelope_smooth = smooth(lever_velocity_envelope,1);

% define active parts of the lever, ignore small movments and
% fill in gaps between likely continuous movements
movethresh = 0.0007;
lever_active = lever_velocity_envelope_smooth > movethresh;

% used to get rid of small movements, but given the new stronger filtering
% anything that crosses threshold is probably usable

% give leeway on both ends to all movements 
movement_leeway = 150; % ms to extend movement total (half on each end)
movement_leeway_filt = ones(movement_leeway,1);
lever_active = logical(conv2(+lever_active,movement_leeway_filt,'same'));

% close gaps of determined size
gap_allowance = 500; %ms

lever_active_switch = diff([0;lever_active;0]);
lever_active_starts = find(lever_active_switch == 1);
lever_active_stops = find(lever_active_switch == -1)-1;

lever_active_movement_times = lever_active_stops - ...
    lever_active_starts;
lever_active_intermovement_times = lever_active_starts(2:end) - ...
    lever_active_stops(1:end-1);

lever_active_fill = lever_active_intermovement_times < gap_allowance;
for i = find(lever_active_fill)'
    lever_active(lever_active_stops(i): ...
        lever_active_starts(i+1)) = 1;
end

% get rid of small movements
% not using for the moment - even small ones seemed real
minimum_movement_fast = 0; %ms (excludes leeway)
minimum_movement = minimum_movement_fast + movement_leeway;

lever_active_switch = diff([0;lever_active;0]);
lever_active_starts = find(lever_active_switch == 1);
lever_active_stops = find(lever_active_switch == -1)-1;

lever_active_movement_times = lever_active_stops - ...
    lever_active_starts;
lever_active_intermovement_times = lever_active_starts(2:end) - ...
    lever_active_stops(1:end-1);

lever_active_erase = lever_active_movement_times < minimum_movement;
for i = find(lever_active_erase)'
    lever_active(lever_active_starts(i): ...
        lever_active_stops(i)) = 0;
end

% the edges of the hilbert envelope always go up, so eliminate first/last
% movements if they're on the edges
lever_active_switch = diff([0;lever_active;0]);
lever_active_starts = find(lever_active_switch == 1);
lever_active_stops = find(lever_active_switch == -1)-1;
if lever_active_starts(1) == 1;
    lever_active(1:lever_active_stops(1)) = 0;
end
if lever_active_stops(end) == length(lever_force_resample)
    lever_active(lever_active_starts(end):end) = 0;
end

%%% MARKED BELOW IS WHAT WAS UPDATED, AP 141027
% Since confident that movement epochs are real, the threshold here is too
% high. It should be the first time it's above 1*standard deviation of
% noise, no high threshold needed.

% refine the lever active starts and stops
lever_startstop = diff([0;lever_active;0]);
lever_active_starts = find(lever_startstop == 1);
lever_active_stops = find(lever_startstop == -1)-1;
noise_cutoff = std(lever_force_resample(~lever_active) - ...
    lever_force_smooth(~lever_active)); % CHANGED AP141027

move_start_values = num2cell(lever_force_smooth(lever_active_starts));
move_start_cutoffs = num2cell(lever_force_smooth(lever_active_starts) ...
    + noise_cutoff);

move_stop_values = num2cell(lever_force_smooth(lever_active_stops));
move_stop_cutoffs = num2cell(lever_force_smooth(lever_active_stops) ...
    + noise_cutoff);

movement_epochs = mat2cell(lever_force_resample(lever_active), ...
    1+lever_active_stops-lever_active_starts,1);

% look for trace consecutively past threshold
thresh_run = 10; % ms CHANGED AP141027

movement_start_offsets = cellfun(@(w,x,y,z) ...
    z:z+find(conv(+(abs(x-w) > abs(y-w)), ...
    ones(thresh_run,1),'same')>=thresh_run,1)- ...
    floor(thresh_run/2),move_start_values,movement_epochs, ...
    move_start_cutoffs,num2cell(lever_active_starts),'uni',false);

movement_stop_offsets = cellfun(@(w,x,y,z) ...
    z-find(conv(+(abs(x(end:-1:1)-w) > abs(y-w)), ...
    ones(thresh_run,1),'same')>=thresh_run,1)+ ...
    floor(thresh_run/2):z,move_stop_values,movement_epochs, ...
    move_stop_cutoffs,num2cell(lever_active_stops),'uni',false);

lever_active(horzcat(movement_start_offsets{:})) = 0;
lever_active(horzcat(movement_stop_offsets{:})) = 0;














