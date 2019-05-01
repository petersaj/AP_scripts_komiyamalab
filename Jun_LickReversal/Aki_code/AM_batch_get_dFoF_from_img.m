function updated = AM_batch_get_dFoF_from_img(mousename,dataroot,n_day1,tmpdir)
%
    updated = false;
            
    if(nargin < 4)
        tmpdir = '/var/tmp/';
    end
    if(nargin < 3)
        n_day1 = 1;
    end
    if(nargin < 2)
        dataroot = fullfile(pwd,'data');
    end
    if(nargin < 1)
        mousename = input('Mouse name?  ','s');
        if(isempty(mousename))
            disp('empty mousename');
            return;
        end
    end
%     
%     mousename = 'JL053';
%     dataroot = fullfile(pwd,'data');

	disp('Listing all ROI files');
    disp(['Data directory is ' dataroot]);

    roi_dir = fullfile(dataroot, mousename, [mousename '_roi']);
    analysis_dir = fullfile(dataroot, mousename, [mousename '_analysis']);
    [s,e]=mkdir(analysis_dir);
    if(s==1 && ~isempty(e))
        disp('Analysis dir already exists, skipping this mouse. Delete the dir to update.');
        return;
    end
    
    roi_dir_names = regexpdir(roi_dir,['.*' mousename '\_\d{6}'],false);
    
    
    if(isempty(roi_dir_names))
        disp('No ROI directories');
        return;
    end
    
%% load ROI files
    disp('Loading ROI files')
    rois = cell(size(roi_dir_names,1),1);
    tif_files = cell(size(roi_dir_names,1),1);
    datestr = cell(size(roi_dir_names,1),1);
    N_day = size(roi_dir_names,1);
    N_loop = zeros(N_day);
    
    tif_file_names = cell(N_day,1);
    roi_file_names = cell(N_day,1);
    
    
    for n_day = n_day1:size(roi_dir_names,1)
        tif_files{n_day} = dir(fullfile(roi_dir_names{n_day},'*.tif'));
        n_loop = 1;
        for n_tif = 1:size(tif_files{n_day},1)
            tif_file_name = tif_files{n_day}(n_tif).name;
            if(~isempty(strfind(tif_file_name,'AVG'))||~isempty(strfind(tif_file_name,'MAX')))
                continue;
            end
            roi_file_name = regexprep(tif_file_name,'\.tif','v2.roi');
            if(~exist(fullfile(roi_dir_names{n_day},roi_file_name),'file'))
                roi_file_name = regexprep(tif_file_name,'\.tif','.roi');
                if(~exist(fullfile(roi_dir_names{n_day},roi_file_name),'file'))
                    disp(['No ROI file for ' tif_file_name]);
                    continue
                end
            end
            roi_file_names{n_day}{n_loop} = roi_file_name;
            tif_file_names{n_day}{n_loop} = tif_file_name;
            rois{n_day}(n_loop) = load(fullfile(roi_dir_names{n_day},roi_file_name),'-mat');
            n_loop = n_loop+1;
        end
        N_loop(n_day) = n_loop-1;
        clear n_loop
        
        
        [fpath, ~, ~] = fileparts(roi_dir_names{n_day});
        datestr{n_day} = regexp(fpath,'\d{6}','match','once');
    end;clear n_day
    
%% calculate dislocation
    disp('Checking dislocation of ROIs');
    for n_day = n_day1:size(roi_dir_names,1)
        
        
        if(isempty(rois{n_day}))
            continue;
        else
            ref_roi = rois{n_day}(1);
        end
        plotted = false(size(ref_roi.polygon.ROI,2),1);
        for n_loop = 1:N_loop(n_day)
            for n_ROI = 1:size(ref_roi.polygon.ROI,2)
                if(length(rois{n_day}(n_loop).polygon.ROI{n_ROI})~=length( ref_roi.polygon.ROI{n_ROI}))
                    disp(['Warning: unmatched number of points to ref, day=' num2str(n_day) ', loop=' num2str(n_loop) ',roi=' num2str(n_ROI)]);
                    if(~plotted(n_ROI))
                        plotted(n_ROI)=true;
                        figure;plot(ref_roi.polygon.ROI{n_ROI}(:,1),ref_roi.polygon.ROI{n_ROI}(:,2),'b' );
                        hold on;
                        plot(rois{n_day}(n_loop).polygon.ROI{n_ROI}(:,1),rois{n_day}(n_loop).polygon.ROI{n_ROI}(:,2),'r' );
                        plot(mean(ref_roi.polygon.ROI{n_ROI}(:,1)),mean(ref_roi.polygon.ROI{n_ROI}(:,2)),'bx' );
                        plot(mean(rois{n_day}(n_loop).polygon.ROI{n_ROI}(:,1)),mean(rois{n_day}(n_loop).polygon.ROI{n_ROI}(:,2)),'r' );
                        title(['unmatched # of points, day=' num2str(n_day) ', loop=' num2str(n_loop) ',roi=' num2str(n_ROI)]);
                        legend({'ref','new'});
                        box off;
                    end
                    rois{n_day}(n_loop).dislocation(n_ROI,:) = mean(rois{n_day}(n_loop).polygon.ROI{n_ROI}) - mean( ref_roi.polygon.ROI{n_ROI});
                else
                    [m, f, ~]=mode(rois{n_day}(n_loop).polygon.ROI{n_ROI} - ref_roi.polygon.ROI{n_ROI});
                    if(f>1/2 * length(ref_roi.polygon.ROI{n_ROI}))
                        rois{n_day}(n_loop).dislocation(n_ROI,:) = m;
                    else
                        disp(['Warning: more than half points are moved, day=' num2str(n_day) ', loop=' num2str(n_loop) ',roi=' num2str(n_ROI)]);
                        rois{n_day}(n_loop).dislocation(n_ROI,:) = mean(rois{n_day}(n_loop).polygon.ROI{n_ROI} - ref_roi.polygon.ROI{n_ROI});
                    end
                end
                
                rois{n_day}(n_loop).dislocation(n_ROI,:) = round(rois{n_day}(n_loop).dislocation(n_ROI,:) );
                
            end
        end
    end

    
     
     
%%
    disp('Loading images');
    
    img_all = cell(size(roi_dir_names,1),1);
    
    for n_day = n_day1:size(roi_dir_names,1)
        disp(['day ' num2str(n_day)]);
        if(isempty(rois{n_day}))
            continue;
        end
        ref_roi = rois{n_day}(1);
        
        %current_position = 0;

        disp('Loading image file info');
        imf = cell(N_loop(n_day),1);
        for n_loop = 1:N_loop(n_day)
            disp(['loop ' num2str(n_loop)]);
            %% get tif file names
%             [~,n,~] = fileparts(tif_file_names{n_day}{n_loop});
%             tif_files{n_day}{n_loop} =  fullfile(roi_dir_names{n_day},[n, '.tif']);
            source = fullfile(roi_dir_names{n_day},tif_file_names{n_day}{n_loop});
            target = fullfile(tmpdir,tif_file_names{n_day}{n_loop});
            copyfile(source,target);
            
            imf{n_loop} = imfinfo(target);
        end
%%        
        n_frames = cellfun(@length,imf);
        
        disp('Loading image file');
        img_all{n_day} = cell(N_loop(n_day),1);
        for n_loop = 1:N_loop(n_day)
            disp(['loop ' num2str(n_loop)]);
            
            for n_frame = 1:size(imf{n_loop},1)
                img = imread(imf{n_loop}(n_frame).Filename,imf{n_loop}(n_frame).Format,'Info',imf{n_loop}(n_frame));
                
                for n_ROI = 1:size(ref_roi.polygon.ROI,2)
                    topleft = max([1 1],floor(min(ref_roi.polygon.bgROI{n_ROI}))-[4 1]);
                    bottomright = min([size(img,2) size(img,1)],ceil(max(ref_roi.polygon.bgROI{n_ROI}))+[4 1]);
                    
                    
                    topleft = topleft+rois{n_day}(n_loop).dislocation(n_ROI,:);
                    bottomright = bottomright+rois{n_day}(n_loop).dislocation(n_ROI,:);
                    
                    for i = 1:2
                        if(topleft(i)<1)
                            bottomright(i) = bottomright(i)-topleft(i)+1;
                            topleft(i) = 1;
                        end
                    end;clear i;
                    for i = 1:2
                        if(bottomright(i)>size(img,3-i))
                            topleft(i) = topleft(i)-bottomright(i)+size(img,3-i) ;
                            bottomright(i) = size(img,3-i);
                        end
                    end;clear i;
                    
                    
                    
                    
                    if(n_loop == 1 && n_frame == 1)
                        img_all{n_day}{n_ROI} = zeros(bottomright(2)-topleft(2)+1,bottomright(1)-topleft(1)+1, sum(n_frames));
                    end
                    
                    img_all{n_day}{n_ROI}(:,:,n_frame+sum(n_frames(1:(n_loop-1)))) = img(topleft(2):bottomright(2),topleft(1):bottomright(1));
                end
                
                if(mod(n_frame,200)==0)
                    disp(['frame ' num2str(n_frame)]);
                end
            end
            %current_position = current_position+size(roi_files{n_day},1);
        end
        
        
        disp('deleting temporary files');
        for n_loop = 1:N_loop(n_day)
            target = fullfile(tmpdir,tif_file_names{n_day}{n_loop});
            delete(target);
        end
        
        disp('Exctracting traces and Saving files');
        roi_bgcorrected_dFoF = zeros(sum(n_frames), size(ref_roi.polygon.ROI,2));
        
        for n_ROI = 1:size(ref_roi.polygon.ROI,2)
            img = img_all{n_day}{n_ROI};
            
            topleft = max([1 1],floor(min(ref_roi.polygon.bgROI{n_ROI}))-[4 1]);
            
            roi = bsxfun(@minus,ref_roi.polygon.ROI{n_ROI},topleft-1);
            bgroi = bsxfun(@minus,ref_roi.polygon.bgROI{n_ROI},topleft-1);
            
            disp(['ROI ' num2str(n_ROI)]);
            [roi_bgcorrected_dFoF(:,n_ROI), intermediates] = AM_get_dFoF_from_img(img, roi, bgroi);
            
            [~,~]=mkdir(fullfile(analysis_dir,[mousename '_' datestr{n_day}]));
            save(fullfile(analysis_dir,[mousename '_' datestr{n_day}], [mousename '_' datestr{n_day} '_' sprintf('%03d',n_ROI) '.mat']),...
                'img','roi','bgroi','intermediates');
        end
        save(fullfile(analysis_dir,[mousename '_' datestr{n_day}], [mousename '_' datestr{n_day} '.mat']),'roi_bgcorrected_dFoF');
        
        
        
%%        
    end
    updated = true;
%%
end