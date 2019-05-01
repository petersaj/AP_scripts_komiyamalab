numframes_diff = cell(size(mice));

for curr_animal = 1:length(mice)
    
    animal = mice(curr_animal).name;

    % Loop through sessions
    for curr_session = 1:length(data_all(curr_animal).im);
                
        curr_concat_data = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});
        
        % Find the summed movie with the selected animal/session
        data_path = ['/usr/local/lab/People/Jun/Data/' animal filesep animal '_roi'];
        data_dir = dir(data_path);
        data_names = {data_dir.name};
        days = cellfun(@(x) x(end-5:end),data_names(3:end),'UniformOutput',false);
        
        curr_abs_session = data_all(curr_animal).session(curr_session);
        
        im_dir = [data_path filesep animal '_' days{curr_abs_session}];
        im_summed_folder = dir([im_dir filesep 'summed*']);
        summed_movie_dir = [im_dir filesep im_summed_folder.name];
        summed_movie_files = dir([summed_movie_dir filesep '*summed*.tif']);
        summed_movie_filename = [summed_movie_dir filesep summed_movie_files.name];
        
        imageinfo = imfinfo(summed_movie_filename,'tiff');
        M = imageinfo(1).Width;
        N = imageinfo(1).Height;
        num_frames = length(imageinfo);
        
        switch mice(curr_animal).scope
            case 1
                sum_frames = 10;
            case 2
                sum_frames = 50;
        end
              
        numframes_diff{curr_animal}(curr_session) = length(curr_concat_data)/sum_frames - ...
            num_frames;
          
    end
    disp(curr_animal)
end


