function AP_fixed_motion_correction_summed_movies(animal,local_copy)
% AP_fixed_motion_correction_summed_movies(animal,local_copy)
%
% After running AP_check_motion_correction, remake summed movies for days
% that had re-corrected frames
% 
% local_copy - if being run on local PC, choose to copy temporary files
% locally instead of reading/writing from server (default TRUE)

%% Set paths 

% Check if local computer
local_comp = ispc;

if nargin == 1 && local_comp
    local_copy = true;
end

if ~local_comp
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    mc_check_filename = ['/usr/local/lab/People/Andy/Data/motion_correction_check' filesep animal '_motion_correction_check.mat'];
else
    data_path = ['Z:\People\Andy\Data\' animal];
    mc_check_filename = ['Z:\People\Andy\Data\motion_correction_check' filesep animal '_motion_correction_check.mat'];
end

data_path_dir = dir(data_path);
data_path_dir_days = cellfun(@(x) any(x),regexp({data_path_dir.name},'^\d\d\d\d\d\d$'));

day_paths = cellfun(@(x) [data_path filesep x],{data_path_dir(([data_path_dir.isdir]) & data_path_dir_days).name},'uni',false);

%% Load the check motion correction file

remake_summed_movie = false(1,length(day_paths));

load(mc_check_filename);
for curr_day = 1:length(day_paths)
    remake_summed_movie(curr_day) = any(~isnan(vertcat( ...
        im_corr(curr_day).postfix{:})));    
end

%% For any days that were re-corrected: re-make summed movie
if any(remake_summed_movie)
    for curr_session = find(remake_summed_movie)
        overwrite = true;    
        AP_downsampleMovie(50,day_paths{curr_session},local_copy,overwrite);
    end
end












