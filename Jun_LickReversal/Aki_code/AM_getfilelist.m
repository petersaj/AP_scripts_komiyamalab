function AM_call_behaviordispatcher(mousename,homedir)
%homedir = pwd;
%mousename = 'JL043';

    dataroot = fullfile(homedir,'data');

    analysis_dir = fullfile(dataroot, mousename, [mousename '_analysis']);
    analysis_file_names = regexpdir(analysis_dir,[mousename '\_\d{6}.*' mousename '\_\d{6}\.mat'],true);

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


    [~, ia, ib] = intersect([analysis_files_all.datenum],[behavior_files_all.datenum]);

    analysis_files = analysis_files_all(ia);
    behavior_files = behavior_files_all(ib);

    [dates, ia, ix] = intersect([analysis_files.datenum],[xsg_dirs_all.datenum]);
    analysis_files = analysis_files(ia);
    xsg_dirs = xsg_dirs_all(ix);

    L = length(dates);
    experiments = struct('datenum',cell(L,1),'day',cell(L,1) ...
        ,'analysis_file',cell(L,1)...
        ,'behavior_file',cell(L,1)...
        ,'xsg_dir',cell(L,1)...
        ,'xsg_files',cell(L,1));

    day0 = min(dates)-1;
    for i = 1:L
        experiments(i).datenum = dates(i);
        experiments(i).day = dates(i) - day0;
        experiments(i).analysis_file = analysis_files(i).fullname;
        experiments(i).behavior_file = behavior_files(i).fullname;
        experiments(i).xsg_dir = xsg_dirs(i).fullname;
        xsg_files = dir(fullfile(xsg_dirs(i).fullname, '*.xsg'));
        experiments(i).xsg_files = arrayfun(@(x)fullfile(experiments(i).xsg_dir,x.name),xsg_files,'UniformOutput',false);
    end
end


    %%
    
    framerate = 6.30040322580645;



    for i = 1:L
        bhv_fullfilename = experiments(i).behavior_file;
        xsg_fullfilenames = experiments(i).xsg_files;
        [experiments(i).bhv_frames, experiments(i).imaged_trials, experiments(i).xsg_trials] ...
            = AP_getBehavior_dispatcher(framerate, bhv_fullfilename, xsg_fullfilenames);        
    end

    

    for i = 1:L
        if(experiments(i).datenum == 735158)
            experiments(i).imaged_trials(experiments(i).xsg_trials == 1)=0;
            experiments(i).xsg_trials(experiments(i).xsg_trials>0) ...
                = experiments(i).xsg_trials(experiments(i).xsg_trials>0)-1;
        end
        
    end
    
    save(mousename, 'experiments');
end

L = length(experiments);
framerate = 6.30040322580645;
%%
for i = 1:L
    varname = [mousename '_' datestr(experiments(i).datenum,'yymmdd') 'long_correct'];
    data = load(experiments(i).analysis_file,varname);
    if(isfield(data,varname))
        disp(['Loading ' experiments(i).analysis_file]);
        matrix = data.(varname);
        experiments(i).used_data = 'long_correct';
    else
        varname = [mousename '_' datestr(experiments(i).datenum,'yymmdd') '_longcorrected_dFoF'];
        data = load(experiments(i).analysis_file,varname);
        if(isfield(data,varname))
            disp('Warning: Data was not found.');
            disp('    Using long_corrected_dFoF instead in ');
            disp(['    ' experiments(i).analysis_file]);
            matrix = data.(varname);
        experiments(i).used_data = 'longcorrected_dFoF';
        else
            disp('Warning: Data was not found.');
            disp('    No long_correct or long_corrected_dFoF found in ');
            disp(['    ' experiments(i).analysis_file]);
            continue;
        end
    end
    
        
    %%
    experiments(i).bhv_frames{1}.states.apply_odor = [0,0];
    experiments(i).imaged_trials(1) = 0;
    warning off;
    bhv_data = load(experiments(i).behavior_file);
    warning on;

    W = length(experiments(i).bhv_frames);  % number of trials

    sec_pre = 3;
    sec_post = 10;

    num_frame_pre = ceil(sec_pre * framerate);
    num_frame_post = ceil(sec_post * framerate);

    experiments(i).trace_pre_odor = NaN(size(matrix,1),num_frame_pre,W);
    experiments(i).trace_post_odor = NaN(size(matrix,1),num_frame_post,W);
    experiments(i).pre_mean = NaN(size(matrix,1),W);
    experiments(i).time_odor_begin = NaN(W,1);

    iti_threshold = 5;

    for w = 1:W
        if(experiments(i).imaged_trials(w))
            experiments(i).time_odor_begin(w) ...
                = experiments(i).bhv_frames{w}.states.apply_odor(1) ...
                + 1000*(experiments(i).xsg_trials(w)-1);
            frame_odor_application = floor(experiments(i).time_odor_begin(w) );
            pre_frames = (frame_odor_application - num_frame_pre):(frame_odor_application-1);
            if(pre_frames(1)<1) 
                disp('Error exceeding matrix size');
                disp(i)
                disp(w)
                disp(pre_frames(1));
                continue;
            end
            post_frames = frame_odor_application:(frame_odor_application + num_frame_post - 1 );
            if(post_frames(end)>size(matrix,2));
                disp('Error exceeding matrix size');
                disp(i)
                disp(w)
                disp(post_frames(end));
                disp(size(matrix,2));
                continue;
            end
            experiments(i).trace_pre_odor(:,:,w) = matrix(:,pre_frames);
            experiments(i).trace_post_odor(:,:,w) = matrix(:,post_frames);
            experiments(i).pre_mean(:,w) = mean(matrix(:,pre_frames),2);
        end
        experiments(i).iti(w) = bhv_data.saved_history.TimesSection_iti{w};
        experiments(i).iti_threshold = iti_threshold;
        experiments(i).iti_ok(w) = bhv_data.saved_history.TimesSection_iti{w} > iti_threshold;

        experiments(i).applied_odor(w) = bhv_data.saved_history.OdorChoiceSection_AppliedOdor{w};

        experiments(i).rewarded_odor(w) = bhv_data.saved_history.OdorChoiceSection_RewardedOdor{w};

        experiments(i).correct_lick(w) = isfield(experiments(i).bhv_frames{w}.states,'correct_lick') ...
            &&  ~isempty(experiments(i).bhv_frames{w}.states.correct_lick);
        experiments(i).incorrect_rejection(w)  ...
            =experiments(i).applied_odor(w) ==  experiments(i).rewarded_odor(w) ...
            & ~experiments(i).correct_lick(w);


        experiments(i).incorrect_lick(w) = isfield(experiments(i).bhv_frames{w}.states,'incorrect_lick') ...
            &&  ~isempty(experiments(i).bhv_frames{w}.states.incorrect_lick);
        experiments(i).correct_rejection(w)  ...
            =experiments(i).applied_odor(w) ~=  experiments(i).rewarded_odor(w) ...
            & ~experiments(i).incorrect_lick(w);


    end
    save([mousename '_trace'], 'experiments');
end

% figure;plot(nanmean(experiments.trace_pre_odor,3)');
% figure;plot(nanmean(experiments.trace_post_odor,3)');



% 'roi_longcorrected_dFoF' is the corrected raw trace I used for my analysis

% code I used to threshold my data, currently I use 3 times STD + mean:
% JL_smooth.m


% This is to visualize all behavior days:
% JunBehaviorPlot.m



% This is to generate data for my bootstrapping approach. DataA is the
% reference data, DataB is covering 6 seconds after odor application
% (roughly 38 frames):
% JL_bootstrappingactivity_6sec.m

% This is the code for bootstrapping:
% HK_bootstrap_onesided.m