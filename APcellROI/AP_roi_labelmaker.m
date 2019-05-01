function varargout = AP_roi_labelmaker(varargin)
% AP_ROI_LABELMAKER MATLAB code for AP_roi_labelmaker.fig
%      AP_ROI_LABELMAKER, by itself, creates a new AP_ROI_LABELMAKER or raises the existing
%      singleton*.
%
%      H = AP_ROI_LABELMAKER returns the handle to a new AP_ROI_LABELMAKER or the handle to
%      the existing singleton*.
%
%      AP_ROI_LABELMAKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AP_ROI_LABELMAKER.M with the given input arguments.
%
%      AP_ROI_LABELMAKER('Property','Value',...) creates a new AP_ROI_LABELMAKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AP_roi_labelmaker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AP_roi_labelmaker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AP_roi_labelmaker

% Last Modified by GUIDE v2.5 05-Nov-2012 17:52:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AP_roi_labelmaker_OpeningFcn, ...
                   'gui_OutputFcn',  @AP_roi_labelmaker_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AP_roi_labelmaker is made visible.
function AP_roi_labelmaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AP_roi_labelmaker (see VARARGIN)

% Choose default command line output for AP_roi_labelmaker
handles.output = hObject;

% Recognize AP_cellROI_pro, if it was opened by that GUI
handles.AP_cellROI_pro = [];
AP_cellROI_pro_input = find(strcmp(varargin, 'AP_cellROI_pro'));
if ~isempty(AP_cellROI_pro_input)
    handles.AP_cellROI_pro = varargin{AP_cellROI_pro_input + 1};
    
    curr_cellroi_handles = guidata(handles.AP_cellROI_pro);
    cellroi_polygon = get(curr_cellroi_handles.roiList,'UserData');
    
    % only populate ROIs if there are currently listed ROIs
    if ~isempty(cellroi_polygon)
        set(handles.roi_list,'string',num2str((1:length(cellroi_polygon.handles))'));
    end   
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AP_roi_labelmaker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AP_roi_labelmaker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in roi_list.
function roi_list_Callback(hObject, eventdata, handles)
% hObject    handle to roi_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roi_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roi_list

listselect = get(handles.roi_list,'value');
selectROI(listselect,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function roi_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in curr_roi_labels.
function curr_roi_labels_Callback(hObject, eventdata, handles)
% hObject    handle to curr_roi_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns curr_roi_labels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from curr_roi_labels


% --- Executes during object creation, after setting all properties.
function curr_roi_labels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curr_roi_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function roi_label_Callback(hObject, eventdata, handles)
% hObject    handle to roi_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi_label as text
%        str2double(get(hObject,'String')) returns contents of roi_label as a double


% --- Executes during object creation, after setting all properties.
function roi_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_label.
function apply_label_Callback(hObject, eventdata, handles)
% hObject    handle to apply_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listselect = get(handles.roi_list,'value');
roi_label = get(handles.roi_label,'string');
num_rois = size(get(handles.roi_list,'string'),1);

if ~isfield(handles,'roi_labels')
    handles.roi_labels = cell(num_rois,1);
end

for i = 1:length(listselect)
    curr_labels = handles.roi_labels{listselect(i)};
    handles.roi_labels{listselect(i)}(length(curr_labels)+1) = {roi_label};
end

% Update handles structure
guidata(hObject, handles);

color_ROI(hObject,eventdata,handles)



% --- Executes on button press in delete_label.
function delete_label_Callback(hObject, eventdata, handles)
% hObject    handle to delete_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

labelselect = get(handles.curr_roi_labels,'value');
listselect = get(handles.roi_list,'value');
num_rois = size(get(handles.roi_list,'string'),1);
applied_labels_all = get(handles.curr_roi_labels,'string');
delete_label = applied_labels_all(labelselect);

for i = 1:length(listselect)
    curr_labels = handles.roi_labels{listselect(i)};
    label_indx = strcmp(curr_labels,delete_label);
    handles.roi_labels{listselect(i)}(label_indx) = [];
end

selectROI(listselect,eventdata,handles)

% Update handles structure
guidata(hObject, handles);

color_ROI(hObject,eventdata,handles)


% --- Executes on button press in load_labels.
function load_labels_Callback(hObject, eventdata, handles)
% hObject    handle to load_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[roilabel_file roilabel_path] = uigetfile('.roilabel','Pick ROIlabel File');
load([roilabel_path roilabel_file],'-MAT')
handles.roi_labels = roi_labels;
set(handles.roi_list,'string',num2str((1:length(roi_labels))'));

listselect = 1;
selectROI(listselect,eventdata,handles)

% Update handles structure
guidata(hObject, handles);

color_ROI(hObject,eventdata,handles)


% --- Executes on button press in save_labels.
function save_labels_Callback(hObject, eventdata, handles)
% hObject    handle to save_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[roilabel_file roilabel_path] = uiputfile('.roilabel','Create save filename');
roi_labels = handles.roi_labels;
save([roilabel_path roilabel_file],'roi_labels')


% --- Executes on button press in load_rois.
function load_rois_Callback(hObject, eventdata, handles)
% hObject    handle to load_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[roi_file roi_path] = uigetfile('.roi','Pick ROI File');
load([roi_path roi_file],'-MAT')
set(handles.roi_list,'string',num2str((1:length(polygon.ROI))'));

% Update handles structure
guidata(hObject, handles);

function selectROI(hObject, eventdata, handles)

listselect = get(handles.roi_list,'value');
if isfield(handles,'roi_labels')
    curr_labels_all = unique(horzcat(handles.roi_labels{listselect}));
    if ~isempty(curr_labels_all)
        set(handles.curr_roi_labels,'string',curr_labels_all,'value',1); 
    else
        set(handles.curr_roi_labels,'string','');
    end
else
    set(handles.curr_roi_labels,'string','');
end


function num_rois_Callback(hObject, eventdata, handles)
% hObject    handle to num_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_rois as text
%        str2double(get(hObject,'String')) returns contents of num_rois as a double

num_rois = str2num(get(handles.num_rois,'string'));
set(handles.roi_list,'string',num2str((1:num_rois)'));

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function num_rois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function select_rois_Callback(hObject, eventdata, handles)
% hObject    handle to select_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of select_rois as text
%        str2double(get(hObject,'String')) returns contents of select_rois as a double


listselect = str2num(get(handles.select_rois,'string'));
if max(listselect) > length(get(handles.roi_list,'string')) || ...
        min(listselect) < 1
    error('Selected ROIs out-of-bounds with number of ROIs')
    return
end
set(handles.roi_list,'value',listselect);
if isfield(handles,'roi_labels')
    curr_labels_all = unique(horzcat(handles.roi_labels{listselect}));
    if ~isempty(curr_labels_all)
        set(handles.curr_roi_labels,'string',curr_labels_all,'value',1); 
    else
        set(handles.curr_roi_labels,'string','');
    end
else
    set(handles.curr_roi_labels,'string','');
end

% --- Executes during object creation, after setting all properties.
function select_rois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function color_ROI(hObject, eventdata, handles)
% find number of current labels, color ROIs in cellROI accordingly
% check for new labels, update label colors
all_labels = [handles.roi_labels{:}]';
all_labels = unique(all_labels);

color_cycle = ['w';'g';'c';'m'];
all_colors = num2cell(color_cycle(rem(0:length(all_labels)-1, ...
    size(color_cycle,1))+1,:));

all_data = [all_labels all_colors];
set(handles.label_colors,'Data',all_data)

if isfield(handles,'AP_cellROI_pro')
    curr_cellroi_handles = guidata(handles.AP_cellROI_pro);
    cellroi_polygon = get(curr_cellroi_handles.roiList,'UserData');
    for curr_label = 1:length(all_labels)
       labelled_rois = ...
           cellfun(@(x) any(strcmp(all_labels{curr_label},x)), ...
           handles.roi_labels);
       set([cellroi_polygon.handles{labelled_rois}],'color', ...
           all_colors{curr_label});
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in label_colors.
function label_colors_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to label_colors (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'AP_cellROI_pro')
    curr_cellroi_handles = guidata(handles.AP_cellROI_pro);
    cellroi_polygon = get(curr_cellroi_handles.roiList,'UserData');
    curr_data = get(handles.label_colors,'Data');
    all_labels = curr_data(:,1);
    all_colors = curr_data(:,2);
    for curr_label = 1:length(all_labels)
       labelled_rois = ...
           cellfun(@(x) any(strcmp(all_labels{curr_label},x)), ...
           handles.roi_labels);
       set([cellroi_polygon.handles{labelled_rois}],'color', ...
           all_colors{curr_label});
    end
end

% Update handles structure
guidata(hObject, handles);


% --- Executes when entered data in editable cell(s) in label_colors.
function label_colors_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to label_colors (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'AP_cellROI_pro')
    curr_cellroi_handles = guidata(handles.AP_cellROI_pro);
    cellroi_polygon = get(curr_cellroi_handles.roiList,'UserData');
    curr_data = get(handles.label_colors,'Data');
    all_labels = curr_data(:,1);
    all_colors = curr_data(:,2);
    for curr_label = 1:length(all_labels)
       labelled_rois = ...
           cellfun(@(x) any(strcmp(all_labels{curr_label},x)), ...
           handles.roi_labels);
       set([cellroi_polygon.handles{labelled_rois}],'Color', ...
           all_colors{curr_label});
    end
end

% Update handles structure
guidata(hObject, handles);
