%% Detect spikes using OOPSI and compare with dot products


animal = 49;
listing_cell = {...
    '120302' ...
    }'

data = struct;
% save figures?
data.save_figures_flag = 1;

for trace_loop_file = 1:length(listing_cell)
    clearvars -except listing_cell trace_loop_file animal data
    curr_loop_folder = ['C:\Users\Andy\Documents\KomiyamaLab\data\' listing_cell{trace_loop_file} '\AP' num2str(animal)];
    % load into structure then unpack to avoid overwriting metaloop vars
    curr_data = [];
    curr_data = load([curr_loop_folder filesep 'Manual_ROI.mat']);
    unstruct(curr_data);
    
    
    
    
    
    %% Vogelstein fast_oopsi
    
    n_best = {};
    P_best = {};
    V = {};
    
    V.dt = 1/30;
    
    for i = 1:size(roi_trace_df,1)
        [n_best{i} P_best{i} V_est]=fast_oopsi(roi_trace_df(i,:),V);
        disp(num2str(i/size(roi_trace_df,1)));
    end
    n_best = [n_best{:}]';
    
    % %% Plot fast_oopsi results
    % roi = 2;
    %
    % figure; hold on;
    % plot(roi_trace_df(roi,:),'k')
    % plot(n_best(roi,:),'r')
    
    %% Correlations via dot product
    % I'm pretty sure this makes sense
    spike_dot = [];
    spike_dot = n_best*n_best';
    self_dot = [];
    self_dot = diag(spike_dot);
    self_dot_rep = [];
    self_dot_rep = repmat(self_dot,1,length(self_dot));
    dot_avg = [];
    dot_avg = (self_dot_rep+self_dot_rep')/2;
    
    spike_integral = [];
    spike_integral = sum(n_best,2)/size(n_best,2);
    % Expected spikes
    spike_expected = [];
    spike_expected = size(n_best,2)*(spike_integral*spike_integral');
    
    norm_spike_dot = [];
    norm_spike_dot = (spike_dot-spike_expected)./dot_avg;
    
    % index and histogram unique dot products
    tril_bin = [];
    tril_bin = ones(size(norm_spike_dot));
    tril_bin = tril(tril_bin,-1);
    tril_bin = logical(tril_bin);
    % figure;
    % hist(norm_spike_dot(tril_bin),20);
    
    data.norm_spike_dot{trace_loop_file} = norm_spike_dot;
    
end