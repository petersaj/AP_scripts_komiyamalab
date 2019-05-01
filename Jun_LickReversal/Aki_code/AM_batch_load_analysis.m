function sessions = AM_batch_load_analysis(sessions)
    % load analysis files
    
    L = size(sessions);
    
    mousename = sessions(1).mousename;
    for i = 1:L
        disp(['Loading ' sessions(i).analysis_file]);
        
        varname = 'roi_bgcorrected_dFoF';
        data = load(sessions(i).analysis_file,varname);
        if(isfield(data,varname))
            sessions(i).matrix = data.(varname);
            sessions(i).used_data = 'roi_bgcorrected_dFoF';
        else
            disp('Warning: Data was not found.');
            disp(['    ' sessions(i).analysis_file]);
        end
%         varname = [mousename '_' datestr(experiments(i).datenum,'yymmdd') 'long_correct'];
%         if(isfield(data,varname))
%             experiments(i).matrix = data.(varname);
%             experiments(i).used_data = 'long_correct';
%         else
%             varname = [mousename '_' datestr(experiments(i).datenum,'yymmdd') '_longcorrected_dFoF'];
%             data = load(experiments(i).analysis_file,varname);
%             if(isfield(data,varname))
%                 disp('Warning: Data was not found.');
%                 disp('    Using longcorrected_dFoF instead in ');
%                 disp(['    ' experiments(i).analysis_file]);
%                 experiments(i).matrix = data.(varname);
%                 experiments(i).used_data = 'longcorrected_dFoF';
%             else
%                 varname = [mousename '_' datestr(experiments(i).datenum,'yymmdd') 'longcorrected_dFoF'];
%                 data = load(experiments(i).analysis_file,varname);
%                 if(isfield(data,varname))
%                     disp('Warning: Data was not found.');
%                     disp('    Using _longcorrected_dFoF instead in ');
%                     disp(['    ' experiments(i).analysis_file]);
%                     experiments(i).matrix = data.(varname);
%                     experiments(i).used_data = 'longcorrected_dFoF';
%                 else
%                     disp('Warning: Data was not found.');
%                     disp('    No long_correct or longcorrected_dFoF or _longcorrected_dFoF found in ');
%                     disp(['    ' experiments(i).analysis_file]);
%                     continue;
%                 end
%             end
%         end
        
        
        

%         sec_pre = 3;
%         sec_post = 10;
% 
%         num_frame_pre = ceil(sec_pre * framerate);
%         num_frame_post = ceil(sec_post * framerate);
%         num_frame = num_frame_pre + num_frame_post;
% 
%         experiments(i).trace= NaN(size(matrix,1),num_frame,W);
%         experiments(i).num_frame_pre = num_frame_pre;
%         experiments(i).num_frame_post = num_frame_post;
%         
%         experiments(i).pre_mean = NaN(size(matrix,1),W);
%         experiments(i).time_odor_begin = NaN(W,1);
%         
%         for w = 1:W
%             if(experiments(i).imaged_trials(w))
%                 experiments(i).time_odor_begin(w) ...
%                     = experiments(i).bhv_frames{w}.states.apply_odor(1) ...
%                     + 1000*(experiments(i).xsg_trials(w)-1);
%                 frame_odor_application = floor(experiments(i).time_odor_begin(w) );
%                 frames = (frame_odor_application - num_frame_pre):(frame_odor_application + num_frame_post - 1 );
%                 if(frames(1)<1 || frames(end)>size(matrix,2)) 
%                     disp('Error exceeding matrix size');
%                     disp(i)
%                     disp(w)
%                     disp(frames);
%                     continue;
%                 end
%                 
%                 experiments(i).trace(:,:,w) = matrix(:,frames);
%                 experiments(i).pre_mean(:,w) = mean(matrix(:,frames(1:num_frame_pre)),2);
%             end
%         end
        
    end
    if(nargout == 0)
        save([mousename '_mat'], 'sessions');
    end
end
    