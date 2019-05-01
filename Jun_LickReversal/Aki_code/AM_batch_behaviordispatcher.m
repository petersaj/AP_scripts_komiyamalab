function sessions = AM_batch_behaviordispatcher(mousename,dataroot)
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
    
   
	disp('Listing all relevant files (analysis, behavior, xsg)');
    disp(['Data directory is ' dataroot]);

    analysis_dir = fullfile(dataroot, mousename, [mousename '_analysis']);
    analysis_file_names = regexpdir(analysis_dir,[mousename '\_\d{6}.*' mousename '\_\d{6}\.mat'],true);
% 
%     analysis_dir = fullfile(dataroot, mousename, [mousename '_roi']);
%     analysis_file_names = regexpdir(analysis_dir,[mousename '\_\d{6}.*analysis\.mat'],true);

    behavior_dir = fullfile(dataroot, mousename, [mousename '_bhv']);
    behavior_file_names = regexpdir(behavior_dir,['.*' mousename '\_\d{6}a\.mat'],true);

    xsg_dir = fullfile(dataroot, mousename, [mousename '_xsg']);
    xsg_dir_names = regexpdir(xsg_dir,['.*' mousename '\_\d{6}'],false);

    %%
    analysis_files_all = struct('fname',cell(length(analysis_file_names),1) ...
                            ,'datestr',cell(length(analysis_file_names),1) ...
                            ,'datenum',cell(length(analysis_file_names),1));
    for i = 1:length(analysis_file_names)
        analysis_files_all(i).fullname = analysis_file_names{i}; 
        
        [~, fname, ~] = fileparts(analysis_files_all(i).fullname);
        analysis_files_all(i).datestr = regexp(fname,'\d{6}','match','once');
%         [fpath, ~, ~] = fileparts(analysis_files_all(i).fullname);
%         analysis_files_all(i).datestr = regexp(fpath,'\d{6}','match','once');
        analysis_files_all(i).datenum = datenum(analysis_files_all(i).datestr,'yymmdd');
    end

    %
    behavior_files_all = struct('fname',cell(length(behavior_file_names),1) ...
                            ,'datestr',cell(length(behavior_file_names),1) ...
                            ,'datenum',cell(length(behavior_file_names),1));
    for i = 1:length(behavior_file_names)
        behavior_files_all(i).fullname = behavior_file_names{i}; 
        [~, fname, ~] = fileparts(behavior_files_all(i).fullname);
        behavior_files_all(i).datestr = regexp(fname,'\d{6}','match','once');
        behavior_files_all(i).datenum = datenum(behavior_files_all(i).datestr,'yymmdd');
    end

    xsg_dirs_all = struct('fname',cell(length(behavior_file_names),1) ...
                            ,'datestr',cell(length(behavior_file_names),1) ...
                            ,'datenum',cell(length(behavior_file_names),1));
    for i = 1:length(xsg_dir_names)
        [xsg_dirs_all(i).fullname,~,~] = fileparts(xsg_dir_names{i}); 
        [~, dname, ~] = fileparts(xsg_dirs_all(i).fullname);
        xsg_dirs_all(i).datestr = regexp(dname,'\d{6}','match','once');
        xsg_dirs_all(i).datenum = datenum(xsg_dirs_all(i).datestr,'yymmdd');
    end

% choose dates which have all necessary files (analysis, behavior, xsg)
    [~, ia, ib] = intersect([analysis_files_all.datenum],[behavior_files_all.datenum]);

    analysis_files = analysis_files_all(ia);
    behavior_files = behavior_files_all(ib);

    [dates, ia, ix] = intersect([analysis_files.datenum],[xsg_dirs_all.datenum]);
    analysis_files = analysis_files(ia);
    xsg_dirs = xsg_dirs_all(ix);
    
%     disp('Analysis files are')
%     disp({analysis_files.fullname})
%     
%     disp('Behavior files are')
%     disp({behavior_files.fullname})
%     
%     disp('XSG directories are')
%     disp({xsg_dirs.fullname})
    
    
    
% store file names
    L = length(dates);
    sessions = struct('datenum',cell(L,1),'day',cell(L,1) ...
        ,'analysis_file',cell(L,1)...
        ,'behavior_file',cell(L,1)...
        ,'xsg_dir',cell(L,1)...
        ,'xsg_files',cell(L,1)...
        ,'mousename',cell(L,1));

    day0 = min(dates)-1;
    for i = 1:L
        sessions(i).datenum = dates(i);
        sessions(i).day = dates(i) - day0;
        sessions(i).analysis_file = analysis_files(i).fullname;
        sessions(i).behavior_file = behavior_files(i).fullname;
        sessions(i).xsg_dir = xsg_dirs(i).fullname;
        xsg_files = dir(fullfile(xsg_dirs(i).fullname, '*.xsg'));
        sessions(i).xsg_files = arrayfun(@(x)fullfile(sessions(i).xsg_dir,x.name),xsg_files,'UniformOutput',false);
        sessions(i).mousename = mousename;
    end
    
    


    %%
    
    framerate = 6.30040322580645;

    disp('Calling AP_getBehavior_dispatcher...')

    for i = 1:L
        disp(['day ' num2str(dates(i) - day0)])
        bhv_fullfilename = sessions(i).behavior_file;
        xsg_fullfilenames = sessions(i).xsg_files;
        [sessions(i).bhv_frames, sessions(i).imaged_trials, sessions(i).xsg_trials] ...
            = AP_getBehavior_dispatcher(framerate, bhv_fullfilename, xsg_fullfilenames);  
        
        % ignore first trial
        sessions(i).bhv_frames{1}.states.apply_odor = [1,1];
        sessions(i).imaged_trials(1) = 0;
        sessions(i).W = length(sessions(i).bhv_frames);  % number of trials
    end

    

%     
% % some specific code for exceptional days
%     for i = 1:L
%         if(strcmp(mousename, 'JL043') && sessions(i).datenum == 735158)
%             disp('Eliminating the first image from day 2 of JL043.')
%             sessions(i).imaged_trials(sessions(i).xsg_trials == 1)=0;
%             sessions(i).xsg_trials(sessions(i).xsg_trials>0) ...
%                 = sessions(i).xsg_trials(sessions(i).xsg_trials>0)-1;
%         end
%         
%     end
%     
%     for i = 1:L
%         if(strcmp(mousename, 'JL049') && sessions(i).datenum == 735257)
%             disp('Eliminating loop 9 from day 1 of JL049.')
%             sessions(i).imaged_trials(sessions(i).xsg_trials == 9)=0;
%             sessions(i).xsg_trials(sessions(i).xsg_trials>9) ...
%                 = sessions(i).xsg_trials(sessions(i).xsg_trials>9)-1;
%         end
%     end

    save([mousename '_dspt'], 'sessions');
    
end