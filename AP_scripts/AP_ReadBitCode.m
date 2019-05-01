function TrialNumber = AP_ReadBitCode(xsg)
% TrialNumber = AP_ReadBitCode(xsg)
%
% 'xsg' input can be either folder name (loads in trial number, assumes
% trace_1 is trial number) or a raw trial number trace
%
% This reads the bitcode from Dispatcher
% This is not to read the default trialnum_indicator with Dispatcher,
% because that doesn't work on our equipment.
% Instead, this reads a bitcode which has a sync signal followed by 12 bits
% for the trial number, which all have 5ms times with 5ms gaps in between.
% Bitcode is most significant bit first (2048 to 1).
% Total time: 5ms sync, 5ms*12 gaps + 5ms*12 bits = 125ms
% The temporal resolution of the linux state machine gives >0.1 ms loss
% per 5 ms period - keep in mind that this causes ~2.6ms to be
% lost over the course of 24 states.
%
% Note: the start of the trial is defined as the START of the bitcode

num_bits = 12;

ThresholdValue = 2;

if ischar(xsg)
    trace = load(xsg,'-mat');
    xsg_sample_rate = trace.header.acquirer.acquirer.sampleRate;
    BinaryThreshold = trace.data.acquirer.trace_1>ThresholdValue;
elseif isvector(xsg)
    xsg_sample_rate = 10000;
    BinaryThreshold = xsg>ThresholdValue;
end

ShiftBinaryThreshold = [NaN; BinaryThreshold(1:end-1)];
% Get raw times for rising edge of signals
rising_bitcode = find(BinaryThreshold==1 & ShiftBinaryThreshold==0);

% Set up the possible bits, 12 values, most significant first
bit_values = [num_bits-1:-1:0];
bit_values = 2.^bit_values;

% Find the sync bitcodes: anything where the difference is larger than the
% length of the bitcode (16 ms - set as 20 ms to be safe)
bitcode_time_samples = 125*(xsg_sample_rate/1000);
bitcode_sync = find(diff(rising_bitcode) > bitcode_time_samples);
% Assume that the first rising edge is a sync signal
if isempty(rising_bitcode)
    
    TrialNumber = [];
    
else
    
    bitcode_sync = rising_bitcode([1;bitcode_sync + 1]);
    % Initialize TrialNumber output - col2 = trial number, col1 = time
    TrialNumber = zeros(length(bitcode_sync),2);
    % for each bitcode sync, check each bit and record as hi or low
    for curr_bitcode_sync = 1:length(bitcode_sync)
        curr_bitcode = zeros(1,num_bits);
        for curr_bit = 1:num_bits
            % boundaries for bits: between the half of each break
            % (bitcode_sync+5ms+2.5ms = 7.5ms)
            bit_boundary_min = bitcode_sync(curr_bitcode_sync) + 7.5*(xsg_sample_rate/1000) + ...
                (curr_bit-1)*10*(xsg_sample_rate/1000);
            bit_boundary_max = bitcode_sync(curr_bitcode_sync) + 7.5*(xsg_sample_rate/1000) + ...
                (curr_bit)*10*(xsg_sample_rate/1000);
            if any(rising_bitcode > bit_boundary_min & rising_bitcode < bit_boundary_max)
                curr_bitcode(curr_bit) = 1;
            end
        end
        curr_bitcode_trial = sum(curr_bitcode.*bit_values);
        % TrialNumber col 2 is trial number
        TrialNumber(curr_bitcode_sync,2) = curr_bitcode_trial;
        % TrialNumber col 1 is time (in seconds)
        TrialNumber(curr_bitcode_sync,1) = bitcode_sync(curr_bitcode_sync)/xsg_sample_rate;
        
        % Catch the rare instance of the xsg file cutting out before the end of the bitcode 
        if bit_boundary_max > length(BinaryThreshold)
           TrialNumber(curr_bitcode_sync,:) = []; 
        end 
    end
    
    % Check here if anything fishy is going on, and warn user
    if ~all(diff(TrialNumber(:,2)))
        warning('TRIAL NUMBER WARNING: Nonconsecutive trials');
    end
    
end