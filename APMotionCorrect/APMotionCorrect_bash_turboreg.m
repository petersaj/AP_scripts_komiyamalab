function APMotionCorrect_bash_turboreg(orig_dir, save_dir)
% Run turboreg from bash

% make save directory if it doesn't exist
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

dir_currfolder = dir(orig_dir);
dir_filenames = {dir_currfolder.name};
tiff_files = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);

tic
master_filelist = dir_filenames(tiff_files);
% append directory
master_filelist = cellfun(@(x) [orig_dir filesep x],master_filelist,'UniformOutput',false);
autoRef = 1;
parsedfile = 1;
channels = 1;
multisession = 0;
create_parsed_avg = 1;

if multisession == 1;
    num_session = max(master_filelist{end}{4});
elseif multisession == 0;
    filelist = master_filelist;
    num_session = 1;
end

for session = 1:num_session
    
    if multisession == 1;
        session_index = cellfun(@(x) x{4}==session,master_filelist);
        filelist = master_filelist(session_index);
    end
    
    % create parsed ref if selected
    if parsedfile == 1 && autoRef == 1;
        APMotionCorrect_parsedRef_bash
    end
    % run turboreg
    APMotionCorrect_turboreg_bash
    %     % get rid of temporary ref files
    %     if exist([savedir filesep 'parsedAverage.tif'],'file')
    %         delete([savedir filesep 'parsedAverage.tif']);
    %     end
    
    % create final average of all files if selected
    if create_parsed_avg == 1
        % get image size from first file
        imageinfo=imfinfo(filelist{1},'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        disp('Creating average from corrected')
        ref_numframes = 0;
        im_avg = zeros(N,M);
        im_max = zeros(N,M);
        % get total number of frames from all files
        for pfiles2 = 1:length(filelist);
            ref_iminfo=imfinfo(filelist{pfiles2},'tiff');
            ref_numframes_long=length(ref_iminfo);
            ref_numframes=ref_numframes + ref_numframes_long/channels;
        end
        for pfiles = 1:length(filelist);
            % Get paths/filenames from GUI (stored in filelist for each
            savedir = save_dir;
            filename = filelist{filenum};
            [filename_path,filename_file,filename_ext] = fileparts(filename);
            % Get saved name of files
            savecorrected = [savedir filesep filename_file '_Turboreg.tif'];
            % Build up average image and max image
            ref_iminfoCurrent=imfinfo(savecorrected,'tiff');
            ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
            for ref_frame = 1:ref_numframesCurrent
                im_build_curr = double(imread(savecorrected,'tiff',ref_frame));
                im_avg(:,:) = im_avg(:,:)+im_build_curr./ref_numframes;
            end
            disp(['Creating average from corrected: ' ...
                num2str(100*(pfiles/length(filelist))) ,'%'])
        end
        disp('Saving average from corrected')
        % save temporary combined average and max
        [filename_path,filename_file,filename_ext] = fileparts(filelist{filenum});
        im_avg = uint16(im_avg);
        for windows_lock = 1:100
            try
                imwrite(im_avg,[savedir filesep filename_file 'AVG.tif'],'tif','Compression','none','WriteMode','overwrite');
                break;
            catch me
                pause(0.2);
            end
            disp('Error! Didn''t Write')
            keyboard
        end
    end
    toc
end