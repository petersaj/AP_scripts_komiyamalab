for filenum = 1:size(filelist,1)
    %% Get initial values, settings, names
    
    % Get paths/filenames from GUI (stored in filelist for each
    savedir = filelist{filenum}{1,1};
    filename = filelist{filenum}{2,1};
    reffilename = filelist{filenum}{3,1};
    [filename_path,filename_file,filename_ext] = fileparts(filename);
    
   
    %% Use appropriate reference
    
    % If turboreg, correct turboreg'd file and use appropriate reference
    if turboreg == 1;
        filename = [savedir filesep filename_file '_Turboreg'];
        channels = 1; % turboreg turns into 1 channel
        if parsedfile == 1 && autoRef == 1;
%             reffilename = [savedir filesep 'parsedTurboregAverage.tif'];
            reffilename = im_t;
        end
    end
    
    % if no turboreg but parsed + autoRef, use parsedAverage reference created earlier
    if parsedfile == 1 && autoRef == 1 && turboreg == 0;
%         reffilename = [savedir filesep 'parsedAverage.tif'];
        reffilename = im_t;
    end
    
    % no parsed, but autoRef, create average image (if turboreg was
    % selected, filename was changed appropriately above)
    if parsedfile == 0 && autoRef == 1;
        imageinfo=imfinfo(filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        im_t = zeros(N,M);
       disp('Creating reference from average image');
        for u = 1:channels:numframes
            im_temp(:,:)=double(imread(filename,'tiff',u));
            im_t = im_t + im_temp./(numframes/channels);
        end
        disp('Done.')
%         imwrite(uint16(im_t),[savedir filesep 'currfileAverage.tif'],'tif','Compression','none','WriteMode','overwrite');
%         reffilename = [savedir filesep 'currfileAverage.tif'];
        reffilename = im_t;
    end
    
           
    %% Run HMM
    
    % make filename for future saved HMM files
    
    if turboreg == 1;
        savecorrected = [savedir filesep filename_file '_Turboreg_HMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' lambda_string];
    else
        savecorrected = [savedir filesep filename_file '_HMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' lambda_string];
    end
    
    % Run HMM
    split_createPI_ManRef(filename,savecorrected,reffilename,0,maxdx,maxdy,channels,lambda);
    load([savecorrected '_PI.mat']);
        
    % Calculate lambda or use given lambda
    % Define the number of frames based on channels
    imageinfo=imfinfo(filename,'tiff');
    numframes=length(imageinfo)/channels;
        
        
%     if autoLambda == 0;
%         lambda = str2num(get(handles.Lambda,'string'));
%         markov_on_PIsave
%     elseif autoLambda == 1;
%         maxlambda;
%     end

    % no one ever used autolambda anyway
    lambda = str2num(get(handles.Lambda,'string'));
    markov_on_PIsave

% Create motion corrected image
    [fixeddata,countdata]=tk_playback_markov2(filename,savecorrected,offsets,edgebuffer,0,lambda,numchannels);
      
    disp(['HMM corrected: ' savecorrected])
    
    if exist([savedir filesep filename_file '.mat'],'file')
        delete([savedir filesep filename_file '.mat']);
    end
    if exist([savedir filesep filename_file '.mat'],'file');
        delete([savedir filesep filename_file '.mat'])
    end
    if turboreg == 1;
        delete([savedir filesep filename_file '_Turboreg.mat'])
        delete([savedir filesep filename_file '_Turboreg.tif'])
    end
    
    % Delete temp PI file
    if exist([savecorrected '_PI.mat'],'file')
        delete([savecorrected '_PI.mat']);
    end
    
    % Save offsets from HMM
    hmm_offsets.y = [];
    hmm_offsets.x = [];
    
    hmm_offsets.y = offsets(1,:);
    hmm_offsets.x = offsets(2,:);
    save([savecorrected '_offsets.mat'],'hmm_offsets');
    
end

% create final average of all files, if picked
create_parsed_avg = get(handles.create_parsed_avg,'value');
if create_parsed_avg == 1;
    % reget channels, may have been changed if Turboreg done before
    channels = get(handles.channels,'value');
    % get image size from first file
    imageinfo=imfinfo(filelist{1}{2,1},'tiff');
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    disp('Creating average from corrected');
    ref_numframes = 0;
    im_avg = zeros(N,M);
    im_max = zeros(N,M);
    % get total number of frames from all files
    for pfiles2 = 1:size(filelist,1);
        ref_iminfo=imfinfo(filelist{pfiles2}{2,1},'tiff');
        ref_numframes_long=length(ref_iminfo);
        ref_numframes=ref_numframes + ref_numframes_long/channels;
    end
    
    for pfiles = 1:size(filelist,1);
        % Get paths/filenames from GUI (stored in filelist for each
        savedir = filelist{filenum}{1,1};
        filename = filelist{filenum}{2,1};
        reffilename = filelist{filenum}{3,1};
        [filename_path,filename_file,filename_ext] = fileparts(filename);
        % Get saved name of files
        if turboreg == 1;
            savecorrected = [savedir filesep filename_file '_Turboreg_HMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' lambda_string '.tif'];
        else
            savecorrected = [savedir filesep filename_file '_HMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' lambda_string '.tif'];
        end
        % Build up average and max image
        ref_iminfoCurrent=imfinfo(filelist{pfiles}{2,1},'tiff');
        ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
        for ref_frame = 1:ref_numframesCurrent
            im_build_curr = double(imread(savecorrected,'tiff',ref_frame));
            im_avg(:,:) = im_avg(:,:)+im_build_curr./ref_numframes;
            im_max_temp(:,:,1) = im_build_curr;
            im_max_temp(:,:,2) = im_max;
            im_max(:,:) = max(im_max_temp,[],3);
        end
        disp(['Creating average from corrected: '...
            num2str(100*(pfiles/size(filelist,1))) ,'%'])
    end
    
    % save temporary combined average
    [filename_path,filename_file,filename_ext] = fileparts(filelist{filenum}{2,1});
    im_avg = uint16(im_avg);
    im_max = uint16(im_max);
    for windows_lock = 1:100
        try
            imwrite(im_avg,[savedir filesep filename_file 'AVG.tif'],'tif','Compression','none','WriteMode','overwrite');
            imwrite(im_max,[savedir filesep filename_file 'MAX.tif'],'tif','Compression','none','WriteMode','overwrite');
            break;
        catch me
            pause(0.2);
        end
    end
end
