% This script is for all types of file manipulation

%% Rename a file
% switch out one specific string in a series of files for another one

old_string = '150209';
new_string = '150219';
filepath = '/usr/local/lab/People/Andy/Data/AP123/150219';

dir_filepath = dir(filepath);
for i = 1:length(dir_filepath)
    str_indx = [];
    str_indx = strfind(dir_filepath(i).name,old_string);
    if ~isempty(str_indx)
       new_filename =  [dir_filepath(i).name(1:str_indx-1) ...
           new_string dir_filepath(i).name(str_indx+length(old_string):end)];
       movefile([filepath filesep dir_filepath(i).name], ...
           [filepath filesep new_filename]);
    end
    disp(num2str(i/length(dir_filepath)))
end

%% Copy over summed movies from server to local

animal = 'AP123';

local_dir = ['C:\Users\Andy\Documents\KomiyamaLab\data\SummedMovies\' animal];
% if this directory doesn't exist, make it
if ~exist(local_dir,'dir')
    mkdir(local_dir);
end

server_dir = ['Z:\People\Andy\Data\' animal];

% look in each folder for a summed movie of the same name
server_dir_dir = dir(server_dir);
% the first two are always . and ..
for i = 3:length(server_dir_dir)
    curr_dir = [server_dir filesep server_dir_dir(i).name filesep 'summed'];
    curr_summed_file = dir([curr_dir filesep '*summed_50.tif']);
    server_sum = [curr_dir filesep ...
        curr_summed_file.name];
    local_sum = [local_dir filesep ...
        curr_summed_file.name];
    % if it's in this folder on the server and not local, copy it over
    if exist(server_sum,'file') && ~exist(local_sum,'file')
        copyfile(server_sum,local_sum);
    end
    disp(['Copied day ' server_dir_dir(i).name])
end
disp('Done')


%% Copy behavior files from ImagingRig3 to motion corrected folders

animal = 'AP164';

raw_path = '/usr/local/lab/Data/ImagingRig3';

data_path = ['/usr/local/lab/People/Andy/Data/' animal];
data_dir = dir(data_path);
data_names = {data_dir.name};
data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
    'UniformOutput',false);
data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
data_folders = data_fullnames(data_isfolders);
data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);

for curr_day = 1:length(days)
    
    day = num2str(days{curr_day});
    curr_source = [raw_path filesep day filesep animal];
    curr_destination = [data_days{curr_day}];
    
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = (dir_filenames(bhv_file == 1));
    
    % in case there's more than one
    for i = 1:length(bhv_filename)
        copyfile([curr_source filesep bhv_filename{i}], ...
            [curr_destination filesep bhv_filename{i}]);
    end
    disp(['Finished day ' num2str(curr_day) '/' num2str(length(days))])

end
%% Copy Simon's behavior files from ImagingRig3 to data folders

animals = {'SC027' 'SC028' 'SC029'};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    raw_path = '/usr/local/lab/Data/ImagingRig3/Simon';
    
    data_path = ['/usr/local/lab/People/Andy/Data/SimonMice_bhv/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    days = cellfun(@(x) x(end-7:end),data_days,'UniformOutput',false);
    
    for curr_day = 1:length(days)
        
        day = num2str(days{curr_day});
        curr_source = [raw_path filesep day filesep];
        curr_destination = [data_days{curr_day}];
        
        dir_currfolder = dir(curr_source);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = (dir_filenames(bhv_file == 1));
        
        % in case there's more than one
        for i = 1:length(bhv_filename)
            if any(strfind(bhv_filename{i},animal))
                copyfile([curr_source filesep bhv_filename{i}], ...
                    [curr_destination filesep bhv_filename{i}]);
            end
        end
        disp(['Finished day ' num2str(curr_day) '/' num2str(length(days))])
        
    end
    
end

%% Copy non-tiff files & subfolders from ImagingRig3 to motion corrected folders

animals = {'AP136'};

for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    raw_path = '/usr/local/lab/Data/ImagingRig3';
    
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);
    
    for curr_day = 1:length(days)
        
        day = num2str(days{curr_day});
        curr_source = [raw_path filesep day filesep animal];
        curr_destination = [data_days{curr_day}];
        
        % Identify non-tiff files / subfolders
        dir_currfolder = dir(curr_source);
        dir_filenames = {dir_currfolder.name};
        nontiff_dir = cellfun(@(x) isempty(strfind(x,'.tif')),dir_filenames);
        nontiff_names = dir_filenames(nontiff_dir);
        
        % copy all the non-tiffs over (first two = ., .., start at 3)
        for i = 3:length(nontiff_names)
            copyfile([curr_source filesep nontiff_names{i}], ...
                [curr_destination filesep nontiff_names{i}]);
        end
        disp(['Finished day ' num2str(curr_day) '/' num2str(length(days))])
        
    end
    
end


%% Copy all XSG/BHV from imaging folder to data folder (non-imaging expts)

animal = 'AP173';
days = {'160109','160110','160111','160112','160113','160114'};

img_path = ['/usr/local/lab/Data/ImagingRig3'];
data_path = ['/usr/local/lab/People/Andy/Data/' animal];

if ~exist(data_path,'dir')
    mkdir(data_path)
end


for curr_day = 1:length(days)
    
    day = days{curr_day};
    curr_source = [img_path filesep day filesep animal];
    curr_destination = [data_path filesep day];
    if ~exist(curr_destination,'dir')
        mkdir(curr_destination)
    end
    
    % bhv files
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = (dir_filenames(bhv_file == 1));
    
    for i = 1:length(bhv_filename)
        copyfile([curr_source filesep bhv_filename{i}], ...
            [curr_destination filesep bhv_filename{i}]);
    end
    
    % xsg files
    dir_currfolder = dir(curr_source);
    dir_filenames = {dir_currfolder.name};
    dir_isfolders = cellfun(@(x) isdir([curr_source filesep x]), ...
        dir_filenames);
    xsg_regexp = regexp(dir_filenames(dir_isfolders),'\w\w\d\d\d\d');
    xsg_folders = cellfun(@(x) ~isempty(x),xsg_regexp);
    dir_folders = dir_filenames(dir_isfolders);
    curr_xsg_folders = dir_folders(xsg_folders);
    
    for i = 1:length(curr_xsg_folders)
        copyfile([curr_source filesep curr_xsg_folders{i}], ...
            [curr_destination filesep curr_xsg_folders{i}]);
    end
    
    
    disp(['Finished day ' num2str(curr_day) '/' num2str(length(days))])

end


%% For Taro Toyoizumi: prepare folders with behavior and dispatcher files
animals = {'AP72' 'AP73' 'AP74' 'AP76' 'AP78' 'AP79'};
for curr_animal = 1:length(animals)
    animal = animals{curr_animal};
    
    % find days
    data_path = ['/usr/local/lab/People/Andy/Data/' animal];
    data_dir = dir(data_path);
    data_names = {data_dir.name};
    data_fullnames = cellfun(@(x) [data_path filesep x],data_names, ...
        'UniformOutput',false);
    data_isfolders = cellfun(@(x) isdir(x),data_fullnames);
    data_folders = data_fullnames(data_isfolders);
    data_regexp = regexp(data_folders,'\d\d\d\d\d\d');
    data_days = data_folders(cellfun(@(x) ~isempty(x),data_regexp));
    days = cellfun(@(x) x(end-5:end),data_days,'UniformOutput',false);
    
    % make folder / copy files for xsg
    mkdir(['/usr/local/lab/People/Andy/Data/' ...
        animal filesep animal '_xsg']);
    for curr_day = 1:length(days)
        mkdir(['/usr/local/lab/People/Andy/Data/' ...
            animal filesep animal '_xsg' filesep ...
            days{curr_day} '_' animal '_xsg']);
        copyfile(['/usr/local/lab/People/Andy/Data/' ...
            animal filesep days{curr_day} filesep 'AP00' animal(3:4)], ...
            ['/usr/local/lab/People/Andy/Data/' ...
            animal filesep animal '_xsg' filesep ...
            days{curr_day} '_' animal '_xsg'])      
    end
    
    % make folder / copy files for dispatcher
    curr_destination = ['/usr/local/lab/People/Andy/Data/' ...
        animal filesep animal '_dispatcher'];
    mkdir(curr_destination);
    for curr_day = 1:length(days)
        curr_source = ['/usr/local/lab/People/Andy/Data/' ...
            animal filesep days{curr_day}];
        dir_currfolder = dir(curr_source);
        dir_filenames = {dir_currfolder.name};
        bhv_file = strncmp('data_@',dir_filenames,6);
        bhv_filename = (dir_filenames(bhv_file == 1));
        
        % in case there's more than one
        for i = 1:length(bhv_filename)
            copyfile([curr_source filesep bhv_filename{i}], ...
                [curr_destination filesep bhv_filename{i}]);
        end
    end
    
    animal
end
































