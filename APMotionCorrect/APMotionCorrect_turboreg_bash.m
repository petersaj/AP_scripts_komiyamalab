if exist('filelist')
    
    % loop through all files in list
    for filenum = 1:length(filelist)
        %% Get paths/filenames from GUI
        
        %savedir = filelist{filenum}{1,1};
        filename = filelist{filenum};
        [filename_path,filename_file,filename_ext] = fileparts(filename);
        
        %filename = get(handles.TIFFFile,'string');
        savedir = save_dir;
        %reffilename = get(handles.TIFFRef,'string');
        
        % Get settings from GUI
        
        savename = [savedir filesep filename_file '_Turboreg'];
        
        imageinfo=imfinfo(filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        %im_s=zeros(N,M,numframes/channels);'
        
        %% Use appropriate reference for autoref
        
        % if parsed + currently autorefd: use avg reference created earlier
        if parsedfile == 1 && autoRef == 1;
%             reffilename = [savedir filesep 'parsedAverage.tif'];
%             im_t = imread(reffilename,'tiff');
            disp('Using parsed average reference')
        % Autoref: make avg image, temp zeros for source image
        elseif parsedfile == 0 && autoRef == 1;
            disp('Creating Reference from Parsed')
            im_t = zeros(N,M);
            for u = 1:channels:numframes
                im_temp(:,:)=double(imread(filename,'tiff',u));
                im_t(:,:) = im_t + im_temp./(numframes/channels);
                disp(['Creating Reference from Parsed: ' ...
                num2str(100*(u/size(numframes,1))) '%']);
            end
            im_s = zeros(N,M);
            
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
        
        Turboreg_AP_pro(im_t,filename,savename,channels);
        
        clear im_s
        
        % Delete extra files, rename finalized one
        disp(['Turboreg Corrected: ' savename])
    end
end

