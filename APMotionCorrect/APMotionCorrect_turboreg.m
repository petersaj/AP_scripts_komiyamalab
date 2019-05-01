if exist('filelist')
    
    % loop through all files in list
    for filenum = 1:size(filelist,1)
        %% Get paths/filenames from GUI
        
        %savedir = filelist{filenum}{1,1};
        filename = filelist{filenum}{2,1};
        reffilename = filelist{filenum}{3,1};
        [filename_path,filename_file,filename_ext] = fileparts(filename);
        
        %filename = get(handles.TIFFFile,'string');
        savedir = get(handles.savedir,'string');
        %reffilename = get(handles.TIFFRef,'string');
        
        % Get settings from GUI
        autoRef = get(handles.AutomaticRef,'value');
        autoLambda = get(handles.AutomaticLambda,'value');
        playback =  get(handles.Playback,'value');
        maxdx = str2num(get(handles.MaxDx,'string'));
        maxdy = str2num(get(handles.MaxDy,'string'));
        channels = get(handles.channels,'value');
        turboreg = get(handles.turboreg,'value');
        lambda = str2num(get(handles.Lambda,'string'));
        parsedfile = get(handles.parsedfile,'value');
        
        savename = [savedir filesep filename_file '_Turboreg'];
        
        imageinfo=imfinfo(filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        %im_s=zeros(N,M,numframes/channels);'
        
        %% Use appropriate reference for autoref
        
        % if parsed + currently autorefd: use avg reference created earlier
        if parsedfile == 1 && autoRef == 1;
            reffilename = [savedir filesep 'parsedAverage.tif'];
%             im_t = imread(reffilename,'tiff');
        
        % Autoref: make avg image, temp zeros for source image
        elseif parsedfile == 0 && autoRef == 1;
            w = waitbar(0,'Turboreg: Creating average image');
            im_t = zeros(N,M);
            for u = 1:channels:numframes
                im_temp(:,:)=double(imread(filename,'tiff',u));
                im_t(:,:) = im_t + im_temp./(numframes/channels);
                waitbar(u/numframes,w);
            end
            im_s = zeros(N,M);
            close(w);
            
        % Get reference or create average for reference
        elseif autoRef == 0;
            im_t = imread(reffilename,'tiff');
        end
        
        % If you've got NaN, make 0 (keep in mind for correlation)
        %im_s(im_s<1)=0;
        
        %% Run Turboreg
        
        % Add ImageJ java controls for turboreg
        spath = javaclasspath('-static');
        spath = cell2mat(spath(1));
        javafolder = strfind(spath,['java' filesep]);
        javafolder = spath(1:javafolder+3);
        
        javaaddpath([javafolder filesep 'jar' filesep 'ij.jar'])
        javaaddpath([javafolder filesep 'jar' filesep 'mij.jar'])
        javaaddpath([javafolder])
        
        Turboreg_AP(im_t,filename,savename,channels,imageinfo);
        
        clear im_s
        
        % Delete extra files, rename finalized one
        disp(['Turboreg Corrected: ' savename])
    end
end

