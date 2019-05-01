function AP_caEvent_vel_threshold_score(data_all,mice,n_animals,n_sessions,n_events)
% AP_caEvent_threshold_score(data_all,mice,n_animals,n_sessions,n_events)
% Manually score calcium events to for "ground truth"
%
% Directions: 
%
% Judge whether "event" shown as red trace is a real calcium event of the
% targeted cell
% 
% Controls:
%
% Mousewheel: scroll through movie
% Left/Right arrows: change event being scored
% Up/Down arrows: score event (green = true, gray = unknown, red = false)
% s: save

%% Prepare data: load images from x sessions from x animals

% Define the time/space to extract for each event
surrounding_time = 10;
cell_border = 100;

total_events = n_animals*n_sessions*n_events;

% Get initials of user (for saving purposes)
user_initials = inputdlg('User initials: ');

% Pick n animals at random (ONLY BSCOPE ANIMALS)
bscope_animals = find([mice.scope] == 2);
use_animals = bscope_animals(randperm(length(bscope_animals),n_animals));

% Initialize the data from each event to extract
events = struct('trace',cell(total_events,1),'std_vel',cell(total_events,1), ...
    'std_df',cell(total_events,1), ...
    'im',cell(total_events,1),'roi',cell(total_events,1), ...
    'maxframe',cell(total_events,1),'activeframes',cell(total_events,1), ...
    'sumframes',cell(total_events,1),'scope',cell(total_events,1), ...
    'animal',cell(total_events,1),'session',cell(total_events,1), ...
    'gad',cell(total_events,1),'eventframes',cell(total_events,1));

% maxframe as [raw max frame, summed max frame]

% Initialize cumulative events index
cumulative_events = 1;
% Loop through animals
disp('Pulling out event information')
for curr_animal_idx = 1:n_animals
    
    curr_animal = use_animals(curr_animal_idx);
    animal = mice(curr_animal).name;
    
    % Pick n sessions at random from current animal
    use_sessions = randperm(length(data_all(curr_animal).im),n_sessions);
    
    % Loop through sessions
    for curr_session_idx = 1:n_sessions
        
        curr_session = use_sessions(curr_session_idx);
        
        % Threshold data, pick n events at random to extract (DON'T USE GAD
        % CELLS AT THE MOMENT)
        curr_concat_data = horzcat(data_all(curr_animal).im(curr_session).roi_trace_df{:});
        [caEventIdxs caEvents caStds_vel caStds_df] = AP_caEvents_vel_threshtest(curr_concat_data,4, ...
            ~data_all(curr_animal).labels.gad);
        event_rois = cellfun(@(x,y) y*ones(size(x)),caStds_vel,num2cell(1:length(caStds_vel))','uni',false);
        event_num = cellfun(@(x) 1:length(x),caStds_vel,'uni',false);
        event_ids = [vertcat(event_rois{:}) horzcat(event_num{:})'];
        use_events = randperm(size(event_ids,1),n_events);
        
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
        
        % Sanity check: if number of summed frames ~= number of frames in
        % trace / sum_frames, don't use this animal (happens at least once)
        trace_sum_frames_check = length(curr_concat_data)/sum_frames == ...
            num_frames;
        if ~trace_sum_frames_check
           disp(['Animal ' animal ' day ' days{curr_abs_session} ...
               ' has unmatching summed movie / raw trace frames, skipping']);
           % delete event slots for the next events
           events(cumulative_events:cumulative_events + n_events - 1) = [];
           continue
        end
        
        % Create TIFF tensor for summed movie
        warning off;
        curr_summed_tensor = TIFFStack(summed_movie_filename);
        warning on;
        
        % Load ROIs
        roi_file = dir([im_dir filesep '*.roi']);
        if length(roi_file) ~= 1
            error('Bad ROI file');
        end
        roi_filename = [im_dir filesep roi_file.name];
        load(roi_filename,'-MAT');
        
        % Loop through chosen events, extract info
        for curr_event_idx = 1:n_events
            
            curr_event = use_events(curr_event_idx);
            curr_event_id = event_ids(curr_event,:);
            
            % Save animal/session/gad
            events(cumulative_events).animal = animal;
            events(cumulative_events).session = curr_abs_session;
            events(cumulative_events).gad = ...
                data_all(curr_animal).labels.gad(curr_event_id(1));
            events(cumulative_events).eventframes = ...
                caEventIdxs{curr_event_id(1)}{curr_event_id(2)};
            
            % Find the frame with the highest df/f, pull out surrounding
            [max_df max_frame_idx] = max(curr_concat_data(caEventIdxs{curr_event_id(1)}{curr_event_id(2)}));
            max_frame = caEventIdxs{curr_event_id(1)}{curr_event_id(2)}(1) + max_frame_idx - 1;
            
            surrounding_frames = ...
                round(data_all(curr_animal).im(curr_session).framerate*surrounding_time);
            extract_frames = max_frame - surrounding_frames : ...
                max_frame + surrounding_frames;
            extract_frames(extract_frames < 1 | extract_frames > length(curr_concat_data)) = [];

            % Save event trace
            events(cumulative_events).trace = curr_concat_data(curr_event_id(1),extract_frames);
            
            % Save which extracted frames were part of this event
            events(cumulative_events).activeframes = ismember(extract_frames, ...
                caEventIdxs{curr_event_id(1)}{curr_event_id(2)});
            
            events(cumulative_events).maxframe(1) = find(extract_frames == max_frame);
            
            % Save scope, get sum frames based on scope
            events(cumulative_events).scope = mice(curr_animal).scope;            
            
            events(cumulative_events).sumframes = sum_frames;
            
            extract_sum_frames = unique(round(extract_frames/sum_frames));
            events(cumulative_events).maxframe(2) = ...
                ceil(events(cumulative_events).maxframe(1)/sum_frames);
          
            % Pull out maximum noise std multiplier for event
            events(cumulative_events).std_vel = caStds_vel{curr_event_id(1)}(curr_event_id(2));
            events(cumulative_events).std_df = caStds_df{curr_event_id(1)}(curr_event_id(2));
            
            % Pull out image of event           
            cell_x = round(cell_border/2);
            cell_y = round(cell_border*(N/M)/2);
            [~,cx,cy] = polycenter(polygon.ROI{curr_event_id(1)}(:,1), ...
                polygon.ROI{curr_event_id(1)}(:,2));
            cx = round(cx);
            cy = round(cy);
            
            cell_y_range = cy-cell_y:cy+cell_y;
            cell_x_range = cx-cell_x:cx+cell_x;
            
            cell_y_use = cell_y_range >= 1 & cell_y_range <= N;
            cell_x_use = cell_x_range >= 1 & cell_x_range <= M;
            
            cell_summed = nan(cell_y*2,cell_x*2,length(extract_sum_frames));
            cell_summed(cell_y_use,cell_x_use,:) = ...
                curr_summed_tensor(cell_y_range(cell_y_use), ...
                cell_x_range(cell_x_use),extract_sum_frames);
            
            % Save event movie
            events(cumulative_events).im = cell_summed;
            
            % Pull out ROI verticies
            events(cumulative_events).roi = polygon.ROI{curr_event_id(1)} - ...
                repmat([cx cy],size(polygon.ROI{curr_event_id(1)},1),1) + ...
                repmat([cell_x+1 cell_y+1],size(polygon.ROI{curr_event_id(1)},1),1);
            
            % Increase running event index
            cumulative_events = cumulative_events + 1;
            
        end       
        disp(['Finished session ' num2str(curr_session_idx) '/' num2str(n_sessions)]);
    end
    disp(['Finished animal ' num2str(curr_animal_idx) '/' num2str(n_animals)]);
end

%% Create GUI for manually scoring activity and compare to thresholds

handles = struct;
% Set the current event
handles.curr_event = 1;

% Create figure
gui_fig = figure;

% Set the axes for plots
handles.im_axes = subplot(7,1,1:6);
handles.trace_axes = subplot(7,1,7);

% Plot the first image
handles.curr_summed_frame = events(1).maxframe(2);
handles.im = imagesc(events(1).im(:,:,handles.curr_summed_frame),'parent',handles.im_axes);
axis(handles.im_axes,'off');
colormap(gray);
set(handles.im_axes,'CLim',[0 max(events(handles.curr_event).im(:))]);

% Plot ROI
handles.roi = line(events(1).roi(:,1),events(1).roi(:,2),'parent',handles.im_axes);

% Plot the first trace
curr_activetrace = events(1).trace;
curr_activetrace(~events(1).activeframes) = NaN;
curr_inactivetrace = events(1).trace;
curr_inactivetrace(events(1).activeframes) = NaN;

hold(handles.trace_axes,'on');
trace_x = linspace(1,size(events(1).im,3),length(events(1).trace));
handles.active_trace_plot = plot(trace_x, ...
    curr_activetrace, ...
    'r','linewidth',1,'parent',handles.trace_axes);
handles.inactive_trace_plot = plot(trace_x, ...
    curr_inactivetrace, ...
    'k','linewidth',1,'parent',handles.trace_axes);

xlim([0 length(events(1).trace)/events(1).sumframes]);
ylim([min(events(1).trace),max(events(1).trace)]);

% Plot lines representing the max frame and the current frame
trace_axes_ylim = get(handles.trace_axes,'ylim');
handles.max_frame_line = line(repmat(events(1).maxframe(2),1,2),trace_axes_ylim, ...
    'color','g','parent',handles.trace_axes,'linewidth',1,'linestyle','--');
handles.curr_frame_line = line(repmat(events(1).maxframe(2),1,2),trace_axes_ylim, ...
    'color','b','parent',handles.trace_axes,'linewidth',1,'linestyle','-');

% Write a title
set(get(handles.im_axes,'Title'),'String','Event 1');

% Package variables for passing to functions
gui_data.handles = handles;
gui_data.events = events;

% Set up scoring
gui_data.score = zeros(size(events));

% Set save name as initials and current time
curr_clock = clock;
gui_data.savename = [user_initials{:} '_' date '_' num2str(ceil(rem(now,1)*10000)) '.eventscores'];

% Store guidata
guidata(gui_fig, gui_data);

% Set functions for mouse wheel / button press and pass necessary handles
set(gui_fig,'KeyPressFcn', {@change_event, gui_fig}); 
set(gui_fig,'WindowScrollWheelFcn',{@slide_movie, gui_fig});

end


%% Slider listener function

function slide_movie(currentObject, eventdata, gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

% Get current event/frame
curr_event = gui_data.handles.curr_event;
old_summed_frame = gui_data.handles.curr_summed_frame;

% Update current frame based on mouse wheel
mouse_wheel_count = eventdata.VerticalScrollCount;
new_summed_frame = old_summed_frame + mouse_wheel_count;
curr_im_size = size(gui_data.events(curr_event).im,3);
if new_summed_frame > curr_im_size
    new_summed_frame = curr_im_size;
elseif new_summed_frame < 1
    new_summed_frame = 1;
end

% Update the image and trace marker
set(gui_data.handles.im,'Cdata',gui_data.events(curr_event).im(:,:,new_summed_frame));
set(gui_data.handles.curr_frame_line,'XData',repmat(new_summed_frame,2,1));

% Update guidata
gui_data.handles.curr_summed_frame = new_summed_frame;
guidata(gui_fig, gui_data);
end


function change_event(currentObject, eventdata, gui_fig)
% Get guidata
gui_data = guidata(gui_fig);

curr_event = gui_data.handles.curr_event;

% Define the keypress events
change_event = false;
score_event = false;
show_roi = false;
save_scores = false;
switch eventdata.Key
    case 'leftarrow'
        if curr_event ~= 1
            new_event = curr_event - 1;
            change_event = true;
            gui_data.handles.curr_event = new_event;
            set(get(gui_data.handles.im_axes,'Title'),'String', ...
                ['Event ' num2str(new_event)]);
        end
    case 'rightarrow'
        if curr_event ~= length(gui_data.events)
            new_event = curr_event + 1;
            change_event = true;
            gui_data.handles.curr_event = new_event;
            set(get(gui_data.handles.im_axes,'Title'),'String', ...
                ['Event ' num2str(new_event)]);
        end
    case 'uparrow'
        if gui_data.score(curr_event) ~= 1
            new_score = gui_data.score(curr_event) + 1;
            gui_data.score(curr_event) = new_score;
            score_event = true;
        end
    case 'downarrow'
        if gui_data.score(curr_event) ~= -1
            new_score = gui_data.score(curr_event) - 1;
            gui_data.score(curr_event) = new_score;
            score_event = true;
        end
    case 'r'
        show_roi = true;
    case 's'
        save_scores = true;
end


% Change events if left/right arrow is pressed
if change_event
    
    % Update background color    
    switch gui_data.score(new_event)
        case -1
            set(gui_fig,'Color',[0.8 0 0])
        case 0 
            set(gui_fig,'Color',[0.8 0.8 0.8]);
        case 1
            set(gui_fig,'Color',[0 0.8 0])
    end
    
    % Update image
    gui_data.handles.curr_summed_frame = gui_data.events(new_event).maxframe(2);
    set(gui_data.handles.im,'Cdata',gui_data.events(new_event).im(:,:,gui_data.handles.curr_summed_frame));
    set(gui_data.handles.im_axes,'CLim',[0 prctile(gui_data.events(new_event).im(:),99)]);
    if all(size(gui_data.events(curr_event).im(:,:,1),1) ~= ...
            size(gui_data.events(new_event).im(:,:,1),1));
        % If new image is different size, replace data and fix x/y limits
        new_im_size = size(gui_data.events(new_event).im(:,:,1));
        set(gui_data.handles.im_axes,'YLim',[0 new_im_size(1)]+0.5, ...
            'XLim',[0 new_im_size(2)]+0.5);
    end
    
    % Update ROI
    set(gui_data.handles.roi,'XData',gui_data.events(new_event).roi(:,1));
    set(gui_data.handles.roi,'YData',gui_data.events(new_event).roi(:,2));
    
    % Update traces 
    curr_activetrace = gui_data.events(new_event).trace;
    curr_activetrace(~gui_data.events(new_event).activeframes) = NaN;
    curr_inactivetrace = gui_data.events(new_event).trace;
    curr_inactivetrace(gui_data.events(new_event).activeframes) = NaN;
    
    trace_x = linspace(1,size(gui_data.events(new_event).im,3),length(gui_data.events(new_event).trace));
    set(gui_data.handles.active_trace_plot,'XData', ...
        trace_x, ...
        'YData',curr_activetrace);
    
    set(gui_data.handles.inactive_trace_plot,'XData', ...
        trace_x, ...
        'YData',curr_inactivetrace);

    set(gui_data.handles.trace_axes,'Xlim',[0 length(gui_data.events(new_event).trace)/gui_data.events(new_event).sumframes]);
    set(gui_data.handles.trace_axes,'Ylim',[min(gui_data.events(new_event).trace),max(gui_data.events(new_event).trace)]);
    
    % Update trace markers
    set(gui_data.handles.curr_frame_line,'XData',repmat(gui_data.handles.curr_summed_frame,2,1));
    set(gui_data.handles.curr_frame_line,'YData',get(gui_data.handles.trace_axes,'YLim'));
    set(gui_data.handles.max_frame_line,'XData',repmat(gui_data.handles.curr_summed_frame,2,1));
    set(gui_data.handles.max_frame_line,'YData',get(gui_data.handles.trace_axes,'YLim'));
    
end


% Score event if up/down arrow is pressed
if score_event
    switch new_score
        case -1
            set(gui_fig,'Color',[0.8 0 0])
        case 0 
            set(gui_fig,'Color',[0.8 0.8 0.8]);
        case 1
            set(gui_fig,'Color',[0 0.8 0])
    end    
end

if show_roi
    % Show/hide blue ROI
    curr_visibility = get(gui_data.handles.roi,'Visible');
    switch curr_visibility
        case 'on'
            new_visibility = 'off';
        case 'off'
            new_visibility = 'on';
    end
    set(gui_data.handles.roi,'Visible',new_visibility);
end


if save_scores
    % Set save directory (common directory for all files)
    save_dir = '/usr/local/lab/People/Andy/Jun_lickreversal/event_detection/vel_threshold/EventScoreDataDrop';
    
    % Set up strucure to save, package data to save
    event_scores = struct('std_vel',cell(size(gui_data.events)),'std_df',cell(size(gui_data.events)), ...
        'animal',cell(size(gui_data.events)), 'session',cell(size(gui_data.events)), ...
        'gad',cell(size(gui_data.events)), 'scope',cell(size(gui_data.events)), ...
        'score',cell(size(gui_data.events)),'eventframes',cell(size(gui_data.events)));
    
    [event_scores.std_vel] = gui_data.events(:).std_vel;
    [event_scores.std_df] = gui_data.events(:).std_df;
    [event_scores.animal] = gui_data.events(:).animal;
    [event_scores.session] = gui_data.events(:).session;
    [event_scores.gad] = gui_data.events(:).gad;
    [event_scores.scope] = gui_data.events(:).scope;
    [event_scores.eventframes] = gui_data.events(:).eventframes;
    
    score_cellarray = num2cell(gui_data.score);
    [event_scores.score] = score_cellarray{:};
    
    % Save the comprehensive structure
    save([save_dir filesep gui_data.savename],'event_scores');
    disp(['Saved: ' gui_data.savename]);
end

% Update guidata
guidata(gui_fig, gui_data);
end
























