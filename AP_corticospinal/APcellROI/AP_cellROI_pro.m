function varargout = AP_cellROI_pro(varargin)
% AP_CELLROI_PRO MATLAB code for AP_cellROI_pro.fig
%      AP_CELLROI_PRO, by itself, creates a new AP_CELLROI_PRO or raises the existing
%      singleton*.
%
%      H = AP_CELLROI_PRO returns the handle to a new AP_CELLROI_PRO or the handle to
%      the existing singleton*.
%
%      AP_CELLROI_PRO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AP_CELLROI_PRO.M with the given input arguments.
%
%      AP_CELLROI_PRO('Property','Value',...) creates a new AP_CELLROI_PRO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AP_cellROI_pro_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AP_cellROI_pro_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AP_cellROI_pro

% Last Modified by GUIDE v2.5 19-Feb-2016 21:52:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AP_cellROI_pro_OpeningFcn, ...
    'gui_OutputFcn',  @AP_cellROI_pro_OutputFcn, ...
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


% --- Executes just before AP_cellROI_pro is made visible.
function AP_cellROI_pro_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AP_cellROI_pro (see VARARGIN)

% Choose default command line output for AP_cellROI_pro
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AP_cellROI_pro wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AP_cellROI_pro_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in createROI.
function createROI_Callback(hObject, eventdata, handles)
% hObject    handle to createROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1) + 1;

% draw ROI polygon
[cellMask,polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2)] = roipoly;

% make the polygon closed if it's not already
if polygon.ROI{n_polygon}(1,1) ~= polygon.ROI{n_polygon}(end,1)
    polygon.ROI{n_polygon} = [polygon.ROI{n_polygon};polygon.ROI{n_polygon}(1,:)];
end

% display just an outline when not selected
polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});


set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);


function tiffFile_Callback(hObject, eventdata, handles)
% hObject    handle to tiffFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tiffFile as text
%        str2double(get(hObject,'String')) returns contents of tiffFile as a double
set(handles.roiList,'String','');
set(handles.roiList,'UserData',[]);

[tiff_fullfilename] = get(handles.tiffFile,'string');
im_s(:,:)=double(imread([tiff_fullfilename],'tiff',1));
axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(im_s,'CDataMapping','scaled');
colormap(gray)
[imgMin imgMax] = caxis;
set(handles.imgMin,'string',imgMin);
set(handles.imgMax,'string',imgMax);
imageinfo=imfinfo(tiff_fullfilename,'tiff');
numframes=length(imageinfo);
set(handles.imgSlider,'Min',1);
if numframes ~= 1;
    set(handles.imgSlider,'Enable','on');
    set(handles.imgSlider,'Max',numframes);
    set(handles.imgSlider,'Value',1);
    set(handles.imgSlider,'SliderStep',[1/(numframes-1), 10/(numframes-1)]);
else
    set(handles.imgSlider,'Enable','off');
end


% --- Executes during object creation, after setting all properties.
function tiffFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tiffFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.roiList,'String','');
set(handles.roiList,'UserData',[]);

% clear old handles if they're around and would cause problems
if isfield(handles,'autosortImg')
    handles = rmfield(handles,'autosortImg');
    handles = rmfield(handles,'guiAxes');
    handles = rmfield(handles,'im');
    guidata(gcbo,handles);
end

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','on');
cd(tiff_path)
if ~iscell(tiff_filename);
    tiff_filename = {tiff_filename};
end

set(handles.tiffFile,'string',[tiff_path tiff_filename{1}]);
% spawn new image window if new or no match
if isfield(handles,'guiAxes') == 0
    % set up figure for keyboard shortcuts
    handles.tiff_fig = figure('KeyPressFcn', {@keyPress, hObject}); 
    % set up the mouse wheel for image scrolling
    set(handles.tiff_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, hObject});
    handles.guiAxes = axes;
    % create scroll bar
    ypos = [0 0 1 0.05];
    handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
    listener_imgSlider = addlistener(handles.imgSlider,'Action',@imgSlider_Listener);
    set(listener_imgSlider,'CallbackTarget',hObject);
    guidata(hObject,handles);
end
if isfield(handles,'guiAxes') == 1
    if any(ismember(findall(0,'type','axes'),handles.guiAxes)) == 0;
        % set up figure for keyboard shortcuts
        handles.tiff_fig = figure('KeyPressFcn', {@keyPress, hObject}); 
        % set up the mouse wheel for image scrolling
        set(handles.tiff_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, hObject});
        handles.guiAxes = axes;
        % create scroll bar
        ypos = [0 0 1 0.05];
        handles.imgSlider = uicontrol('style','slider','units','normalized','position',ypos,'min',0,'max',1,'value',0);
        listener_imgSlider = addlistener(handles.imgSlider,'Action',@imgSlider_Listener);
        set(listener_imgSlider,'CallbackTarget',hObject);
        guidata(hObject,handles);
    end
end

% Get number of channels and channel to load
channels = str2num(get(handles.channels,'string'));

% initialize concatenated matrix by finding total number of frames
total_numframes = 0;
for i = 1:length(tiff_filename);
    img_filename = [tiff_path tiff_filename{i}];
    imageinfo = [];
    imageinfo=imfinfo([img_filename],'tiff');
    numframes=length(imageinfo);
    total_numframes = total_numframes + numframes;
end

total_numframes = floor(total_numframes/channels(1));

M=imageinfo(1).Width;
N=imageinfo(1).Height;
im_concat = zeros(N*M,total_numframes,'uint16');

for i = 1:length(tiff_filename);
    
    img_filename = [tiff_path tiff_filename{i}];
    
    %%%% load in image file
    imageinfo = [];
    imageinfo=imfinfo([img_filename],'tiff');
    numframes=length(imageinfo);
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    im = [];
    im_temp = [];
    disp('Loading file....')
    loadframes = channels(2):channels:numframes;
        for loadframe_idx = 1:length(loadframes);
            loadframe = loadframes(loadframe_idx);
            im_temp = imread([img_filename],'tiff',loadframe,'Info',imageinfo);
            im_temp = im_temp(:);
            im_concat(:,loadframe_idx) = im_temp;
            disp(['Frame ' num2str(loadframe_idx) '/' num2str(total_numframes)]);
        end
end

handles.im = im_concat;
% disp('CURRENTLY LOADING IN AS DOUBLE')
% handles.im = double(handles.im);
%%%%

axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
xlim([0 M]);
ylim([0 N]);
colormap(gray)
%[imgMin imgMax] = caxis;
% make 5% of the pixels the highest

caxis([0 2000]);
set(handles.imgMin,'string',0);
set(handles.imgMax,'string',2000);
numframes=size(handles.im,2);
set(handles.imgSlider,'Min',1);
if numframes ~= 1;
    set(handles.imgSlider,'Enable','on');
    set(handles.imgSlider,'Max',numframes);
    set(handles.imgSlider,'Value',1);
    set(handles.imgSlider,'SliderStep',[1/(numframes-1), 10/(numframes-1)]);
else
    set(handles.imgSlider,'Enable','off');
end

% store the image size
handles.N = N;
handles.M = M;

% make sure figure toolbar is available
if isfield(handles,'tiff_fig')
    set(handles.tiff_fig,'toolbar','figure');
end

% Update handles structure
guidata(hObject, handles);


function imgMin_Callback(hObject, eventdata, handles)
% hObject    handle to imgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgMin as text
%        str2double(get(hObject,'String')) returns contents of imgMin as a double

imgMin = get(handles.imgMin,'string');
imgMax = get(handles.imgMax,'string');
imgMin = str2num(imgMin);
imgMax = str2num(imgMax);
axes(handles.guiAxes);
caxis([imgMin imgMax]);


% --- Executes during object creation, after setting all properties.
function imgMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function imgMax_Callback(hObject, eventdata, handles)
% hObject    handle to imgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgMax as text
%        str2double(get(hObject,'String')) returns contents of imgMax as a double

imgMin = get(handles.imgMin,'string');
imgMax = get(handles.imgMax,'string');
imgMin = str2num(imgMin);
imgMax = str2num(imgMax);
axes(handles.guiAxes);
caxis([imgMin imgMax]);

% --- Executes during object creation, after setting all properties.
function imgMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in roiList.
function roiList_Callback(hObject, eventdata, handles)
% hObject    handle to roiList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns roiList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roiList
axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
listselect = get(handles.roiList,'value');

selectROI(polygon.handles{listselect},eventdata,handles)

% --- Executes during object creation, after setting all properties.
function roiList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roiList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in copyROI.
function copyROI_Callback(hObject, eventdata, handles)
% hObject    handle to copyROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.guiAxes);
listselect = get(handles.roiList,'value');
polygon = get(handles.roiList,'UserData');
polygon.ROI{listselect} = getPosition(polygon.handles{listselect});
n_polygon = size(get(handles.roiList,'string'),1) + 1;

polygon.ROI{n_polygon} = polygon.ROI{listselect};
polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1)+2,polygon.ROI{n_polygon}(:,2)+2);

set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);

% record old position, turn into line
oldROI = listselect;
polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
delete(polygon.handles{oldROI})
polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});

% select new ROI
set(handles.roiList,'value',n_polygon);
set(handles.roiList,'UserData',polygon);
selectROI(polygon.handles{n_polygon},eventdata,handles)


% --- Executes on button press in getTrace.
function getTrace_Callback(hObject, eventdata, handles)
% hObject    handle to getTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listselect = get(handles.roiList,'value');
polygon = get(handles.roiList,'UserData');
try
    polygon.ROI{listselect} = getPosition(polygon.handles{listselect});
catch me
end
AP_cellROI_getTrace_pro
axes(handles.cellTrace);
hold off;
plot(cellTrace);
xlabel('Frame');
ylabel('{\Delta}F/F')


% --- Executes on button press in deleteROI.
function deleteROI_Callback(hObject, eventdata, handles)
% hObject    handle to deleteROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.guiAxes);
listselect = get(handles.roiList,'value');
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1);

delete(polygon.handles{listselect});
polygon.ROI(listselect) = [];
polygon.handles(listselect) = [];

set(handles.roiList,'value',1);
set(handles.roiList,'string',num2str((1:n_polygon-1)'));
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in avgImg.
function avgImg_Callback(hObject, eventdata, handles)
% hObject    handle to avgImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of avgImg
avgImg = get(handles.avgImg,'value');

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

if avgImg == 0;
    set(handles.imgSlider,'Enable','On');
    tiff_fullfilename = get(handles.tiffFile,'string');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
    
elseif avgImg == 1;
    set(handles.imgSlider,'Enable','Off');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(nanmean(handles.im,2),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on slider movement.
function imgSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imgSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% % get frame number from slider, round appropriately
% tiff_fullfilename = get(handles.tiffFile,'string');
% frame = get(handles.imgSlider,'Value');
% frame = round(frame);
% set(handles.imgSlider,'Value',frame);
%
% [tiff_path tiff_filename] = fileparts(tiff_fullfilename);
% imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
% M=imageinfo(1).Width;
% N=imageinfo(1).Height;
%
% % get and plot corresponding image
% set(handles.imgFrame,'String',num2str(frame));
%
% axes(handles.guiAxes);
% imgMin = str2num(get(handles.imgMin,'string'));
% imgMax = str2num(get(handles.imgMax,'string'));
% im_h = findobj('type','image');
% delete(im_h(1));
% hold on;
% imagesc(reshape(handles.im(:,frame),N,M),'CDataMapping','scaled');
% colormap(gray);
% caxis = ([imgMin imgMax]);
% set(gca,'Children',flipud(get(gca,'Children')));
%
% % if behavior plot is open, update current line
% showBehavior = get(handles.showBehavior,'value');
% if showBehavior == 1;
%     axes(handles.behaviorPlot);
%     hold on;
%     if isfield(handles,'behaviorLine')
%         delete(handles.behaviorLine);
%     end
%     handles.behaviorLine = line([frame frame],[1 4],'color','red');
%     xlim([frame-100 frame+100]);
%     % update handles
%     guidata(hObject, handles);
% end

% % redraw polygons
% axes(handles.guiAxes);
% polygon = get(handles.roiList,'UserData');
% n_polygon = size(get(handles.roiList,'string'),1);
% for redraw_polygon = 1:n_polygon;
%     polygon.ROI{redraw_polygon} = getPosition(polygon.handles{redraw_polygon});
%     impoly(gca,[polygon.ROI{redraw_polygon}(:,1) polygon.ROI{redraw_polygon}(:,2)]);
%     setColor(polygon.handles{redraw_polygon},'g')
% end



% --- Executes during object creation, after setting all properties.
function imgSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% give this slider real-time/continuous callback
listener_imgSlider = addlistener(hObject,'Action',@imgSlider_Listener);


function imgSlider_Listener(hObject, eventdata, handles)

handles = guidata(hObject);

% get frame number from slider, round appropriately
tiff_fullfilename = get(handles.tiffFile,'string');
frame = get(handles.imgSlider,'Value');
frame = round(frame);
set(handles.imgSlider,'Value',frame);
%
% [tiff_path tiff_filename] = fileparts(tiff_fullfilename);
% imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
% M=imageinfo(1).Width;
% N=imageinfo(1).Height;
M = handles.M;
N = handles.N;
numframes = size(handles.im,2);

% get and plot corresponding image
set(handles.imgFrame,'String',num2str(frame));

axes(handles.guiAxes);
imgMin = str2num(get(handles.imgMin,'string'));
imgMax = str2num(get(handles.imgMax,'string'));
im_h = findobj('type','image');
delete(im_h(1));
hold on;
imagesc(reshape(handles.im(:,frame),N,M),'CDataMapping','scaled');
colormap(gray);
caxis = ([imgMin imgMax]);
set(gca,'Children',flipud(get(gca,'Children')));

% if behavior plot is open, update current line
showBehavior = get(handles.showBehavior,'value');
if showBehavior == 1;
    axes(handles.behaviorPlot);
    hold on;
    if isfield(handles,'behaviorLine')
        delete(handles.behaviorLine);
    end
    curr_frame = frame*(length(handles.lever_force_resample)/numframes);
    handles.behaviorLine = line([curr_frame curr_frame],ylim,'color','k');
    set(handles.behaviorPlot, ...
        'XLim',[curr_frame-200 curr_frame+200]);
    set(handles.behaviorPlot, 'YLim', [min(handles.lever_force_resample) ...
        max(handles.lever_force_resample) ]);
    % update handles
    guidata(hObject, handles);
end

axes(handles.guiAxes);

function imgSlider_MouseWheel(currentObject, eventdata, hObject)
% executes when mouse wheel is scrolled in figure
handles = guidata(hObject);

wheel_clicks = eventdata.VerticalScrollCount;

numframes = size(handles.im,2);

% get frame number from slider, round appropriately
tiff_fullfilename = get(handles.tiffFile,'string');
frame = get(handles.imgSlider,'Value');
frame = round(frame)+wheel_clicks;

% constrain scroll movement to number of frames
if frame < 1 || frame > numframes
    return
end

set(handles.imgSlider,'Value',frame);

M = handles.M;
N = handles.N;

% get and plot corresponding image
set(handles.imgFrame,'String',num2str(frame));

axes(handles.guiAxes);
imgMin = str2num(get(handles.imgMin,'string'));
imgMax = str2num(get(handles.imgMax,'string'));
im_h = findobj('type','image');
delete(im_h(1));
hold on;
imagesc(reshape(handles.im(:,frame),N,M),'CDataMapping','scaled');
colormap(gray);
caxis = ([imgMin imgMax]);
set(gca,'Children',flipud(get(gca,'Children')));

% if behavior plot is open, update current line
showBehavior = get(handles.showBehavior,'value');
if showBehavior == 1;
    axes(handles.behaviorPlot);
    hold on;
    if isfield(handles,'behaviorLine')
        delete(handles.behaviorLine);
    end
    handles.behaviorLine = line([frame frame],[1 4],'color','red');
    xlim([frame-100 frame+100]);
    % update handles
    guidata(hObject, handles);
end

function imgFrame_Callback(hObject, eventdata, handles)
% hObject    handle to imgFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgFrame as text
%        str2double(get(hObject,'String')) returns contents of imgFrame as a double

frame = str2num(get(handles.imgFrame,'String'));
set(handles.imgSlider,'Value',frame);

% --- Executes during object creation, after setting all properties.
function imgFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autosortON.
function autosortON_Callback(hObject, eventdata, handles)
% hObject    handle to autosortON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
autosortON = get(handles.autosortON,'value');
polygon = get(handles.roiList,'UserData');

if autosortON == 1;
    %if isfield(polygon,'autosort') && ~isfield(handles,'autosortImg');
    if isfield(handles,'autosortImg');
        axes(handles.autosortImg);
        cla(handles.autosortImg);
        linkaxes([handles.guiAxes handles.autosortImg],'xy');
        guidata(gcbo,handles);
        imagesc(sum(polygon.autosort,3));
        colormap(gray);
    else
        axes(handles.guiAxes);
        handles.autosortImg = axes('Position',get(handles.guiAxes,'Position'))
        linkaxes([handles.guiAxes handles.autosortImg],'xy');
        guidata(gcbo,handles);
        imagesc(sum(polygon.autosort,3));
        colormap(gray);
    end
    %     elseif isfield(polygon,'autosort') && isfield(handles,'autosortImg');
    %         axes(handles.autosortImg);
    %         if isempty(findobj('type','image'))
    %             imagesc(sum(polygon.autosort,3));
    %         end
    %     end
else
    if isfield(handles,'autosortImg');
        axes(handles.guiAxes);
    end
end


% Hint: get(hObject,'Value') returns toggle state of autosortON



function smwidth_Callback(hObject, eventdata, handles)
% hObject    handle to smwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smwidth as text
%        str2double(get(hObject,'String')) returns contents of smwidth as a double


% --- Executes during object creation, after setting all properties.
function smwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double


% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function arealims_Callback(hObject, eventdata, handles)
% hObject    handle to arealims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of arealims as text
%        str2double(get(hObject,'String')) returns contents of arealims as a double


% --- Executes during object creation, after setting all properties.
function arealims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to arealims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nPCs_Callback(hObject, eventdata, handles)
% hObject    handle to nPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nPCs as text
%        str2double(get(hObject,'String')) returns contents of nPCs as a double


% --- Executes during object creation, after setting all properties.
function nPCs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxrounds_Callback(hObject, eventdata, handles)
% hObject    handle to maxrounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxrounds as text
%        str2double(get(hObject,'String')) returns contents of maxrounds as a double


% --- Executes during object creation, after setting all properties.
function maxrounds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxrounds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function termtol_Callback(hObject, eventdata, handles)
% hObject    handle to termtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of termtol as text
%        str2double(get(hObject,'String')) returns contents of termtol as a double


% --- Executes during object creation, after setting all properties.
function termtol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to termtol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mu_Callback(hObject, eventdata, handles)
% hObject    handle to mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mu as text
%        str2double(get(hObject,'String')) returns contents of mu as a double


% --- Executes during object creation, after setting all properties.
function mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autosortGO.
function autosortGO_Callback(hObject, eventdata, handles)
% hObject    handle to autosortGO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tiff_fullfilename = get(handles.tiffFile,'string');
nPCs = str2num(get(handles.nPCs,'string'));
mu = str2num(get(handles.mu,'string'));
termtol = str2num(get(handles.termtol,'string'));
maxrounds = str2num(get(handles.maxrounds,'string'));
dt =str2num(get(handles.dt,'string'));
smwidth = str2num(get(handles.smwidth,'string'));
thresh = str2num(get(handles.thresh,'string'));
% arealims: 1 num gives minimum
arealims = str2num(get(handles.arealims,'string'));


% select all files to perform this on
[tiff_filename,tiff_path]=uigetfile('*.tif','Choose TIFF files','Multiselect','on');
if ~iscell(tiff_filename)
    tiff_filename = {tiff_filename};
end
% store ROIs in roilist edit box
h_autosort = figure;
AP_cellROIAutosort_pro;
close(h_autosort);
ica_segments = permute(ica_segments,[2 3 1]);
ica_results.segments = ica_segments;
ica_results.sig = ica_sig;
set(handles.roiList,'UserData',ica_results);

% % create polygons from ica segments (not trivial, fix later maybe?)
% for ica_seg = 1:size(ica_results.segments,3)
%     ica_seg_edge = edge(ica_results.segments(:,:,ica_seg));
%     ica_seg_edge_indx = find(ica_seg_edge == 1);
%     [y x] = ind2sub([N M],ica_seg_edge_indx);
%     polygon.handles{ica_seg} = impoly(gca,[x y]);
%     setColor(polygon.handles{ica_seg},'g')
% end

% draw edges from ica segments
axes(handles.guiAxes)
for ica_seg = 1:size(ica_results.segments,3)
    tiff_fullfilename = get(handles.tiffFile,'string');
    imageinfo=imfinfo([tiff_fullfilename],'tiff');
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    first_nonzero = find(ica_results.segments(:,:,ica_seg) > 0,1);
    [y_nonzero x_nonzero] = ind2sub([N M],first_nonzero);
    ica_edge = bwtraceboundary(ica_results.segments(:,:,ica_seg),[y_nonzero x_nonzero],'N');
    polygon.handles{ica_seg} = line(ica_edge(:,2),ica_edge(:,1));
    set(polygon.handles{ica_seg}, 'ButtonDownFcn', {@selectROI, handles});
    % get x and y right for later
    polygon.ROI{ica_seg} = [ica_edge(:,2),ica_edge(:,1)];
    hold on;
end
polygon.autosort = ica_results.segments;
set(handles.roiList,'string',num2str([1:size(ica_results.segments,3)]'));
set(handles.roiList,'UserData',polygon);



% --- Executes on button press in saveROI.
function saveROI_Callback(hObject, eventdata, handles)
% hObject    handle to saveROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of saveROI
tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[roi_savefile roi_savepath] = uiputfile('.roi','Create save filename',tiff_filename);

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

polygon = get(handles.roiList,'UserData');

% make sure positions are updated in case one is currently selected
axes(handles.guiAxes)
oldpoly = cellfun('isclass', polygon.handles,'impoly');
if any(oldpoly)
    oldROI = find(oldpoly == 1);
    polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
    delete(polygon.handles{oldROI})
    polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
    set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
end

% update polygon handles
set(handles.roiList,'UserData',polygon);

% don't save the trace
%AP_cellROI_saveTrace

% don't bother saving handles
polygon = rmfield(polygon,'handles');
if isfield(polygon,'bghandles');
    polygon = rmfield(polygon,'bghandles')
end

save([roi_savepath roi_savefile],'polygon')%'roi_trace',


% --- Executes on button press in loadROI.
function loadROI_Callback(hObject, eventdata, handles)
% hObject    handle to loadROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[roi_file roi_path] = uigetfile('.roi','Pick ROI File',tiff_filename);

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

axes(handles.guiAxes);

% delete existing polygons and lines
currpoly = findobj(handles.guiAxes,'tag','impoly');
delete(currpoly)
currpoly = findobj(handles.guiAxes,'type','line');
delete(currpoly)

load([roi_path roi_file],'-MAT')

for n_polygon = 1:length(polygon.ROI)
    
    % make the polygon closed if it's not already
    if polygon.ROI{n_polygon}(1,1) ~= polygon.ROI{n_polygon}(end,1)
        polygon.ROI{n_polygon} = [polygon.ROI{n_polygon};polygon.ROI{n_polygon}(1,:)];
    end
    
    polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2));
    set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
end

set(handles.roiList,'Value',1);
set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);


function selectROI(hObject, eventdata, handles)
% Execute when clicking on ROI polygon

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');

oldpoly = cellfun('isclass', polygon.handles,'impoly');
if any(oldpoly)
    oldROI = find(oldpoly == 1);
    polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
    % make the polygon closed if it's not already
    if polygon.ROI{oldROI}(1,1) ~= polygon.ROI{oldROI}(end,1)
        polygon.ROI{oldROI} = [polygon.ROI{oldROI};polygon.ROI{oldROI}(1,:)];
    end
    delete(polygon.handles{oldROI})
    polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
    set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
end


% redraw all polygons in blue
currpoly = findobj(handles.guiAxes,'type','line');

for n_polygon = 1:length(polygon.ROI)
    set(currpoly(n_polygon), 'color', 'b');
end

% draw selected polygon in red
roiClicked = find([polygon.handles{:}] == hObject);
delete(polygon.handles{roiClicked})
if length(polygon.ROI{roiClicked}) < 20
    polygon.handles{roiClicked} = ...
        impoly(gca,[polygon.ROI{roiClicked}(:,1) polygon.ROI{roiClicked}(:,2)]);
    setColor(polygon.handles{roiClicked},'r')
else
    polygon.handles{roiClicked} = line(polygon.ROI{roiClicked}(:,1), polygon.ROI{roiClicked}(:,2),'color','r');
    set(polygon.handles{roiClicked}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'value',roiClicked);
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in roi_right.
function roi_right_Callback(hObject, eventdata, handles)
% hObject    handle to roi_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');

for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    polygon.ROI{i}(:,1) = polygon.ROI{i}(:,1) + 1+(x10*10);
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = 1+(x10*10);
    dy = 0;
    theta = 0;
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end

% --- Executes on button press in roi_left.
function roi_left_Callback(hObject, eventdata, handles)
% hObject    handle to roi_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');

for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    polygon.ROI{i}(:,1) = polygon.ROI{i}(:,1) - 1-(x10*10);
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = -1-(x10*10);
    dy = 0;
    theta = 0;
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end

% --- Executes on button press in roi_up.
function roi_up_Callback(hObject, eventdata, handles)
% hObject    handle to roi_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');

for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    polygon.ROI{i}(:,2) = polygon.ROI{i}(:,2) - 1-(x10*10);
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = 0;
    dy = -1-(x10*10);
    theta = 0;
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end

% --- Executes on button press in roi_down.
function roi_down_Callback(hObject, eventdata, handles)
% hObject    handle to roi_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');

for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    polygon.ROI{i}(:,2) = polygon.ROI{i}(:,2) + 1+(x10*10);
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = 0;
    dy = 1+(x10*10);
    theta = 0;
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end



% --- Executes on button press in roi_rotate_r.
function roi_rotate_r_Callback(hObject, eventdata, handles)
% hObject    handle to roi_rotate_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
aspect_ratio = M/N;

axes(handles.guiAxes);
half_x = max(floor(xlim./2));
half_y = max(floor(ylim./2));
center = [half_x half_x 0];
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');
rotate_matrix = makehgtform('zrotate',(pi/180)+(x10*10*(pi/180)));
translate_matrix = makehgtform('translate',center);
% scale because of 4x columns than rows
scale_factor = aspect_ratio;

for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    for curr_point = 1:length(polygon.ROI{i})
        old_point = [polygon.ROI{i}(curr_point,1) ...
            polygon.ROI{i}(curr_point,2)*aspect_ratio 0 1];
        new_point = translate_matrix*rotate_matrix...
            *inv(translate_matrix)*old_point';
        new_point = [new_point(1) new_point(2)/aspect_ratio];
        polygon.ROI{i}(curr_point,:) = new_point;
    end
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = 0;
    dy = 0;
    theta = (pi/180)+(x10*10*(pi/180));
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end

% --- Executes on button press in roi_rotate_l.
function roi_rotate_l_Callback(hObject, eventdata, handles)
% hObject    handle to roi_rotate_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
aspect_ratio = M/N;

axes(handles.guiAxes);
half_x = max(floor(xlim./2));
half_y = max(floor(ylim./2));
center = [half_x half_x 0];
polygon = get(handles.roiList,'UserData');
x10 = get(handles.roi_x10,'Value');
% create transform matricies
rotate_matrix = makehgtform('zrotate',(-pi/180)-(x10*10*(pi/180)));
translate_matrix = makehgtform('translate',center);
% scale because of 4x columns than rows
scale_factor = aspect_ratio;

% transform polygons point by point, redraw and update
for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    for curr_point = 1:length(polygon.ROI{i})
        old_point = [polygon.ROI{i}(curr_point,1) ...
            polygon.ROI{i}(curr_point,2)*aspect_ratio 0 1];
        new_point = translate_matrix*rotate_matrix...
            *inv(translate_matrix)*old_point';
        new_point = [new_point(1) new_point(2)/aspect_ratio];
        polygon.ROI{i}(curr_point,:) = new_point;
    end
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

% if aligning by image, do the same to the image
if ~isempty(findobj('type','figure','name','im_align'));
    dx = 0;
    dy = 0;
    theta = -(pi/180)-(x10*10*(pi/180));
    
    % center image and constrain output to image size
    [N M] = size(handles.im_align);
    
    udata = [1 N] - (N+1)/2;% input image x
    vdata = [1 M] - (M+1)/2;% input image y
    
    xdata = udata;% output image x
    ydata = vdata;% output image y
    
    % make transform matrix
    A = ...
        [cos(theta) sin(theta) 0;
        -sin(theta) cos(theta) 0;
        dx          dy          1];
    
    tform_rotate = maketform('affine',A);
    
    im_align_tform = imtransform(handles.im_align,tform_rotate,...
        'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    
    axes(handles.align);
    clf
    handles.align = axes;
    imagesc(handles.align_orig);
    hold on;
    handles.im_align = im_align_tform;
    handles.im_align_h = imagesc(handles.im_align);
    set(handles.im_align_h,'AlphaData',0.5);
    guidata(hObject, handles);
end

% --- Executes on button press in roi_x10.
function roi_x10_Callback(hObject, eventdata, handles)
% hObject    handle to roi_x10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_x10

function keyPress(currentObject, eventdata, hObject)
% use shortcuts when image is selected
handles = guidata(hObject);
switch eventdata.Key
    case 'r'
        createROI_Callback(hObject,eventdata,handles);
    case 'h'
        helpROI_Callback(hObject, eventdata, handles);
    case 'delete'
        deleteROI_Callback(hObject, eventdata, handles);
end



% --- Executes on button press in roi_auto_align.
function roi_auto_align_Callback(hObject, eventdata, handles)
% hObject    handle to roi_auto_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

polygon = get(handles.roiList,'UserData');

align_image = get(handles.align_image,'value');
fine_correct = get(handles.fine_correct,'value');

[tiff_fullfilename] = get(handles.tiffFile,'string');
imageinfo=imfinfo(tiff_fullfilename,'tiff');
polygon = get(handles.roiList,'UserData');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

aspect_ratio = M/N;

im_curr_range1 = input('First usable frame from today? (#) ');
im_curr_range2 = input('Last usable frame from today? (#/end) ','s');

if strmatch(im_curr_range2,'end') == 1
   im_curr_range2 = numframes;
else
   im_curr_range2 = str2num(im_curr_range2);
end

if im_curr_range2 > numframes
    im_curr_range2 = numframes;
end

im_curr_range = im_curr_range1:im_curr_range2;
im_curr_avg = reshape(mean(handles.im(:,im_curr_range),2),N,M,1);

[orig_file orig_path] = uigetfile('*.tif','Choose original image');
imageinfo_orig = imfinfo([orig_path orig_file]);
im_orig_avg = double(zeros(N,M));
disp('Creating average from reference day...')
for curr_frame = 1:length(imageinfo_orig)
    im_orig = double(imread([orig_path orig_file],'tiff',curr_frame));
    im_orig_avg = im_orig_avg + im_orig/length(imageinfo_orig);
end
disp('Aligning ROIs...')

% Estimate average ROI radius
polygon_areas = zeros(length(polygon.ROI),1);
for curr_roi = 1:length(polygon.ROI)
    polygon_areas(curr_roi) = polyarea( ...
        polygon.ROI{curr_roi}(:,1),aspect_ratio*polygon.ROI{curr_roi}(:,2));
end

cell_radius = round(sqrt(mean(polygon_areas)/pi));

% Luminance normalize average images
h = fspecial('average', [cell_radius ceil(cell_radius/aspect_ratio)]);
im_orig_avg_lumnorm = im_orig_avg./imfilter(im_orig_avg,h);
im_curr_avg_lumnorm = im_curr_avg./imfilter(im_curr_avg,h);

% Occasionally get NaN/Inf around edges, set to zero
im_orig_avg_lumnorm(isinf(im_orig_avg_lumnorm) | ...
    isnan(im_orig_avg_lumnorm)) = 0;
im_curr_avg_lumnorm(isinf(im_curr_avg_lumnorm) | ...
    isnan(im_curr_avg_lumnorm)) = 0;

% Go through each ROI, figure out the best new center from the old center
% this assumes little rotation and no movement of rois from straight load
s = 5;
s_aspect = round(s/aspect_ratio);

aligned_polygon = polygon.ROI;
uncorrected_flag = false(length(aligned_polygon),1);
corr_offset_all = zeros(length(polygon.ROI),2);

% estimate centers of all ROIs
roi_centers = zeros(length(polygon.ROI),2);
for curr_roi = 1:length(polygon.ROI)
    % estimate center of roi
    center_x = round(mean(polygon.ROI{curr_roi}(:,1)));
    center_y = round(mean(polygon.ROI{curr_roi}(:,2)));
    roi_centers(curr_roi,:) = [center_x center_y];
end

for curr_roi = 1:length(polygon.ROI)
    % estimate center of roi
    center_x = roi_centers(curr_roi,1);
    center_y = roi_centers(curr_roi,2);
    roi_snip = [];
    if center_y-s_aspect*cell_radius <= 0 || center_x-s*cell_radius <= 0 || ...
            center_y+s_aspect*cell_radius > N || center_x+s*cell_radius > M
        uncorrected_flag(curr_roi) = true;
        continue
    end
    % pull out s*cell_radius around the cell, cross correlate
    roi_snip = im_orig_avg_lumnorm(center_y-s_aspect*cell_radius:center_y+s_aspect*cell_radius, ...
        center_x-s*cell_radius:center_x+s*cell_radius);
    cc = normxcorr2(roi_snip,im_curr_avg_lumnorm);
    
    % get the point of max correlation
    xcorr_midpoint = round(length(cc)/2);
    [max_c,imax] = max(cc(:));
    [ypeak xpeak] = ind2sub(size(cc),imax(1));
    corr_offset = [(xpeak-s*cell_radius)-center_x ...
        (ypeak-s_aspect*cell_radius)-center_y];
    
    corr_offset_all(curr_roi,:) = corr_offset;
    
    % record the new positions based on center offset
    aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
        repmat(corr_offset,size(aligned_polygon{curr_roi},1),1);
end

% find abnormal offsets, reset them and flag them as uncorrected
corr_offset_all_dist = sqrt(corr_offset_all(:,1).^2 + ...
    (s_aspect*corr_offset_all(:,2)).^2);

% use 4x the distance of spread of middle 80% as offset cutoff
offset_prctile = prctile(corr_offset_all_dist(corr_offset_all_dist ~= 0),[10 90]);
offset_cutoff = 4*diff(offset_prctile);

abnormal_offsets = find(corr_offset_all_dist < ...
    median(corr_offset_all_dist) - offset_cutoff | ...
    corr_offset_all_dist > ...
    median(corr_offset_all_dist) + offset_cutoff);
aligned_polygon(abnormal_offsets) = polygon.ROI(abnormal_offsets);
uncorrected_flag(abnormal_offsets) = true;

% need at least 100 frames for this
if fine_correct == 1 && length(im_curr_range) >= 100;
    % Fine correct based on activity if possible
    disp('Fine aligning by activity...')
    fine_limit = cell_radius;
    finecorrected_flag = false(length(aligned_polygon),1);
    % use the luminance normalized one for this
    im_lum = imfilter(im_orig_avg,h);
    im_lum = im_lum(:);
    for curr_roi = 1:length(polygon.ROI)
        % get polygon mask, strip down to just cell
        temp_mask = poly2mask(aligned_polygon{curr_roi}(:,1)', ...
            aligned_polygon{curr_roi}(:,2)',N,M);
        curryx = [floor(min(aligned_polygon{curr_roi}(:,2))) ...
            floor(min(aligned_polygon{curr_roi}(:,1))); ...
            ceil(max(aligned_polygon{curr_roi}(:,2))) ...
            ceil(max(aligned_polygon{curr_roi}(:,1)))];
        % if the roi goes outside the frame, drop it
        if any(curryx(1,:) <= 0) || any(curryx(2,:) > [N M]);
            continue
        end
        
        % Apparently it doesn't matter whether it's luminance normalized or
        % not. I haven't double checked that this makes mathematical sense.
        
        % get the correlation between the values in the area
        % use luminance normalized
        %     curr_trace = mean(double(handles.im(temp_mask(:),im_curr_range))./ ...
        %         repmat(im_lum(temp_mask(:)),1,length(im_curr_range)));
        % for non-luminance normalized:
        curr_trace = mean(double(handles.im(temp_mask(:),im_curr_range)));
        curryx_fine = curryx + [-fine_limit -fine_limit; fine_limit fine_limit];
        % if the search is outside the frame, drop it
        if any(curryx_fine(1,:) <= 0) || any(curryx_fine(2,:) > [N M]);
            continue
        end
        curr_cc = zeros(N,M);
        for curr_y = curryx_fine(1):curryx_fine(2)
            for curr_x = curryx_fine(3):curryx_fine(4)
                curr_ind = sub2ind([N M],curr_y,curr_x);
                % do the correlation with a luminance normalized trace
                %             temp_cc = corrcoef(curr_trace, ...
                %                 double(handles.im(curr_ind,im_curr_range))./ ...
                %                 repmat(im_lum(curr_ind),1,length(im_curr_range)));
                temp_cc = corrcoef(curr_trace, ...
                    double(handles.im(curr_ind,im_curr_range)));
                curr_cc(curr_y,curr_x) = temp_cc(2);
            end
        end
        
        % threshold the correlation by what's in the roi range
        norm_roi_thresh = mean(curr_cc(temp_mask(:))) ...
            -1*std(curr_cc(temp_mask(:)));
        
        % if the threshold is less than zero, drop it
        if norm_roi_thresh <=0
            continue
        end
        
        norm_roi = im2bw(curr_cc,norm_roi_thresh);

        % if nothing there, drop it
        if sum(norm_roi(:)) == 0
            continue
        end
        
        % if the size is significantly different (+/- 25%), drop it
        % Also, find center of gravity for previously aligned ROI:
        % using filexchange script 'polycenter' (BiomeCardio)
        [area aligned_center_x aligned_center_y] = ...
            polycenter(aligned_polygon{curr_roi}(:,1), ...
            aligned_polygon{curr_roi}(:,2));
        
        if sum(norm_roi(:)) > area*1.25 || sum(norm_roi(:)) < area*0.75
            continue
        end
        
        % find center of gravity for over threshold points
        [thresh_y thresh_x] = ind2sub([N M],find(norm_roi));
        thresh_center_y = mean(thresh_y);
        thresh_center_x = mean(thresh_x);
        
        % find offset based on the two centers aligning
        fine_offset_y = round(thresh_center_y - aligned_center_y);
        fine_offset_x = round(thresh_center_x - aligned_center_x);
        fine_offset = [fine_offset_x fine_offset_y];
        
        % realign ROIs
        aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
            repmat(fine_offset,size(aligned_polygon{curr_roi},1),1);
        
        finecorrected_flag(curr_roi) = true;
    end
end

% THIS DIDN'T WORK: Try instead getting the ROI shape via corrcoef, or just
% move the center based on the center middle of the corrcoef region

% % TEST: try fine-tuning alignment based on centering over active region
% % try doing this via std: build up
% disp('Fine-correcting ROIs by activity')
% im_std_lumnorm = zeros(N,M);
% for i = 1:length(im_curr_range)
%     curr_im_lumnorm = (imfilter(reshape(mean(handles.im(:,i),2),N,M,1),h) ...
%         ./im_curr_avg - im_curr_avg_lumnorm);   
%     im_std_lumnorm = im_std_lumnorm+curr_im_lumnorm;   
% end
% 
% % the move shouldn't be more than half of a cell
% fine_limit = cell_radius;
% finecorrected_flag = false(length(aligned_polygon),1);
% for curr_roi = 1:length(polygon.ROI)
%     % get polygon mask, strip down to just cell
%     temp_mask = poly2mask(aligned_polygon{curr_roi}(:,1)', ...
%         aligned_polygon{curr_roi}(:,2)',N,M);
%     curryx = [floor(min(aligned_polygon{curr_roi}(:,2))) ...
%         floor(min(aligned_polygon{curr_roi}(:,1))); ...
%         ceil(max(aligned_polygon{curr_roi}(:,2))) ...
%         ceil(max(aligned_polygon{curr_roi}(:,1)))];
%     % if the roi goes outside the frame, drop it
%     if any(curryx(1,:) < 0) || any(curryx(2,:) > [N M]);
%         continue
%     end
%     temp_mask = temp_mask(curryx(1):curryx(2),curryx(3):curryx(4));
%     curr_conv = conv2(im_std_lumnorm,double(temp_mask),'same');
%     curryx_fine = curryx + [-fine_limit -fine_limit; fine_limit fine_limit];
%     % if the search is outside the frame, drop it
%     if any(curryx_fine(1,:) < 0) || any(curryx_fine(2,:) > [N M]);
%         continue
%     end
%     curr_conv_limit = zeros(size(curr_conv));
%     curr_conv_limit(curryx_fine(1):curryx_fine(2), ...
%         curryx_fine(3):curryx_fine(4)) = ...
%         curr_conv(curryx_fine(1):curryx_fine(2), ...
%         curryx_fine(3):curryx_fine(4));
%     [std_y std_x] = ind2sub(size(curr_conv_limit), ...
%         find(curr_conv_limit == max(curr_conv_limit(:)),1));
%     curr_mask_center = round(mean(curryx));
%     curr_activity_offset = [std_y std_x] - curr_mask_center;
%     % if it moves more than the fine limit, drop it
%     if any(sqrt(sum(curr_activity_offset.^2)) > fine_limit)
%         continue
%     end
%     finecorrected_flag(curr_roi) = true;
%     aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
%         repmat(curr_activity_offset,size(aligned_polygon{curr_roi},1),1);
% end

% Estimate offset for unaligned ROIs
median_corr_offset = median(corr_offset_all(~uncorrected_flag));
for curr_roi = find(uncorrected_flag)'
   % estimate offset based on nearby ROIs, get distance
   curr_roi_center = [];
   curr_roi_center = repmat(roi_centers(curr_roi,:),length(roi_centers),1);
   curr_roi_dist = [];
   curr_roi_dist = sqrt(sum((roi_centers - curr_roi_center).^2,2));
   curr_roi_dist = [curr_roi_dist (1:length(curr_roi_dist))'];
   % don't use ROIs that weren't aligned
   curr_roi_dist = curr_roi_dist(~uncorrected_flag,:);
   % sort them by distance, use the top 3
   curr_roi_dist = sortrows(curr_roi_dist,1);
   
   corr_offset = median(corr_offset_all(curr_roi_dist(1:3,2),:));
   
   corr_offset_all(curr_roi,:) = corr_offset;
   
   aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
        repmat(corr_offset,size(aligned_polygon{curr_roi},1),1);
end

figure
title('Aligned ROIs')
hold on
ylim([0 N]);
xlim([0 M]);
set(gca,'YDir','reverse')
set(gca,'color','k');
% plot the old ROIs in dark gray
for curr_roi = 1:length(polygon.ROI)
    roi1(curr_roi) = patch(polygon.ROI{curr_roi}(:,1),polygon.ROI{curr_roi}(:,2),[0.3 0.3 0.3]);
end

% plot the new ROIs in light gray
for curr_roi = 1:length(aligned_polygon)
    roi2(curr_roi) = patch(aligned_polygon{curr_roi}(:,1),aligned_polygon{curr_roi}(:,2),[0.7 0.7 0.7]);
end

% set the uncorrected ROIs to red
set(roi2(uncorrected_flag),'FaceColor','red')

% set the fine corrected ROIs to green
if exist('finecorrected_flag','var')
    set(roi2(finecorrected_flag),'FaceColor','green')
end

% draw lines between old and new ROIs, label ROI numbers
for curr_roi = 1:length(polygon.ROI)
    center_x1 = round(mean(polygon.ROI{curr_roi}(:,1)));
    center_y1 = round(mean(polygon.ROI{curr_roi}(:,2)));
    
    center_x2 = round(mean(aligned_polygon{curr_roi}(:,1)));
    center_y2 = round(mean(aligned_polygon{curr_roi}(:,2)));
    
    line([center_x1 center_x2],[center_y1 center_y2], ...
        'color','w','linewidth',2)
    roi_numbers_h(curr_roi) = text(center_x1,center_y1,num2str(curr_roi), ...
        'color','white','FontWeight','bold','FontSize',14);
end
disp('Done.')
apply_rois = input('Apply these changes? (y/n) ','s');
if strmatch(apply_rois,'y') == 1
    % save over the old positions
    polygon.ROI = aligned_polygon;
    axes(handles.guiAxes);

    % delete existing polygons and lines
    currpoly = findobj(handles.guiAxes,'tag','impoly');
    delete(currpoly)
    currpoly = findobj(handles.guiAxes,'type','line');
    delete(currpoly)
    
    % redraw the lines in their new positions
    for n_polygon = 1:length(polygon.ROI)
        polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2));
        set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
    end
    
    % update roi list
    set(handles.roiList,'UserData',polygon);
end



% --- Executes on button press in roi_numbers.
function roi_numbers_Callback(hObject, eventdata, handles)
% hObject    handle to roi_numbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of roi_numbers

roi_numbers = get(handles.roi_numbers,'value');
polygon = get(handles.roiList,'UserData');
axes(handles.guiAxes);

if roi_numbers == 1;
    roi_centers = zeros(length(polygon.ROI));
    for i = 1:length(polygon.ROI);
        % get center x and y from polygon boundaries
        center_x = mean(polygon.ROI{i}(:,1));
        center_y = mean(polygon.ROI{i}(:,2));
        roi_centers(i,1) = center_x;
        roi_centers(i,2) = center_y;
    end
    for i = 1:length(polygon.ROI);
        roi_numbers_h(i) = text(roi_centers(i,1),roi_centers(i,2),num2str(i),'color','red');
    end
    handles.roi_numbers_h = roi_numbers_h;
    guidata(hObject, handles);
elseif roi_numbers == 0;
    delete(handles.roi_numbers_h)
    handles = rmfield(handles,'roi_numbers_h');
    guidata(hObject, handles);
end


% --- Executes on button press in show_roi.
function show_roi_Callback(hObject, eventdata, handles)
% hObject    handle to show_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_roi
axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
show_roi = get(handles.show_roi,'value');

if show_roi == 1;
    oldpoly = cellfun('isclass', polygon.handles,'impoly');
    if any(oldpoly)
        oldROI = find(oldpoly == 1);
        polygon.ROI{oldROI} = getPosition(polygon.handles{oldROI});
        delete(polygon.handles{oldROI})
        polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
        set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
    end
    for currpoly = 1:length(polygon.handles);        
        set(polygon.handles{currpoly},'Visible','on');
        if isfield(handles,'roi_numbers_h')
            try
            set(handles.roi_numbers_h(currpoly),'Visible','on');
            end
        end
    end
else
    for currpoly = 1:length(polygon.handles);
        set(polygon.handles{currpoly},'Visible','off');
        if isfield(handles,'roi_numbers_h')
            try
            set(handles.roi_numbers_h(currpoly),'Visible','off');
            end
        end
    end
end

set(handles.roiList,'UserData',polygon);


% --- Executes on button press in align_image.
function align_image_Callback(hObject, eventdata, handles)
% hObject    handle to align_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
align_image = get(handles.align_image,'value');
[tiff_fullfilename] = get(handles.tiffFile,'string');
imageinfo=imfinfo(tiff_fullfilename,'tiff');
polygon = get(handles.roiList,'UserData');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

[orig_file orig_path] = uigetfile('*.tif','Choose original image');
imageinfo_orig = imfinfo([orig_path orig_file]);
im_orig_avg = double(zeros(N,M));
disp('Creating average from reference day...')
for curr_frame = 1:length(imageinfo_orig)
    im_orig = double(imread([orig_path orig_file],'tiff',curr_frame));
    im_orig_avg = im_orig_avg + im_orig/length(imageinfo_orig);
end
disp('Done.')

% scale images to be 0-255 uint 8
% im_orig = uint8((im_orig./(max(im_orig(:))))*255);
% im_align = reshape(mean(handles.im,2),N,M);
% im_align = uint8((im_align./(max(im_align(:))))*255);
% assume ~2000 is the max
im_orig_avg = uint8((im_orig_avg./(2000)*255));
im_align = reshape(mean(handles.im,2),N,M);
im_align = uint8((im_align./(2000)*255));

[input_points base_points] = cpselect(im_orig_avg,im_align,'Wait',true);
h = cp2tform(input_points,base_points,'nonreflective similarity');

% account for rows ~= columns
aspect_ratio = M/N;

axes(handles.guiAxes);
% transform polygons point by point, redraw and update
for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    for curr_point = 1:length(polygon.ROI{i})
        old_point = [polygon.ROI{i}(curr_point,1) ...
            polygon.ROI{i}(curr_point,2)*aspect_ratio 0 1];
        new_point = tformfwd(h,old_point(1),old_point(2));
        new_point = [new_point(1) new_point(2)/aspect_ratio];
        polygon.ROI{i}(curr_point,:) = new_point;
    end
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);

guidata(hObject, handles);


% --- Executes on button press in stdImg.
function stdImg_Callback(hObject, eventdata, handles)
% hObject    handle to stdImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stdImg

stdImg = get(handles.stdImg,'value');

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

if stdImg == 0;
    set(handles.imgSlider,'Enable','On');
    tiff_fullfilename = get(handles.tiffFile,'string');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
    colormap(gray)
    caxis([0 2000]);
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
    
elseif stdImg == 1;
    set(handles.imgSlider,'Enable','Off');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    
    % Build std up frame by frame (faster, more memory efficient)
    im_avg = mean(handles.im,2);    
    im_std = zeros(size(im_avg));
    for i = 1:length(size(handles.im,2))
        curr_im_diff = (double(handles.im(:,i)) - im_avg).^2;
        im_std = im_std+curr_im_diff./size(handles.im,2);
    end
    im_std = sqrt(im_std);

    imagesc(reshape(im_std,N,M),'CDataMapping','scaled');
    colormap(gray)
    imgMin = min(im_std(:));
    imgMax = max(im_std(:))/2;
    caxis([imgMin imgMax]);
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on button press in helpROI.
function helpROI_Callback(hObject, eventdata, handles)
% hObject    handle to helpROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This attempts to find the shape of the cell given the center and the
% approximate size of a cell using correlation 

cell_diameter = str2num(get(handles.cell_diameter,'string'));
% make cell diameter an even number for ease in halving
cell_diameter = 2*(ceil(cell_diameter/2));

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

aspect_ratio = M/N;

cell_radius_x = cell_diameter/2;
cell_radius_y = round(cell_radius_x/aspect_ratio);
cell_area = round(pi*cell_radius_x*cell_radius_y);

axes(handles.guiAxes);
[all_x all_y] = ginput(1);
% reinitialize the mouse wheel, some weird matlab java glitch
set(handles.tiff_fig,'WindowScrollWheelFcn',{@imgSlider_MouseWheel, hObject});
x = round(all_x);
y = round(all_y);
center_indx = sub2ind([N M],y,x);
corrcoef_grid = [];

[cell_pixels_x cell_pixels_y] = meshgrid...
    (x-cell_radius_x:x+cell_radius_x, ...
    y-cell_radius_y:y+cell_radius_y);

% get rid of any pixels out of bounds
cell_pixels_x = cell_pixels_x(:);
cell_pixels_y = cell_pixels_y(:);
oob = false(size(cell_pixels_x));
oob = cell_pixels_x <= 0 | cell_pixels_y <= 0 |...
    cell_pixels_x > M | cell_pixels_y > N;
cell_pixels_x(oob) = [];
cell_pixels_y(oob) = [];

cell_estimate_indx = sub2ind([N M],cell_pixels_y,cell_pixels_x);
curr_trace = mean(handles.im(cell_estimate_indx,:));

search_factor = 4;
[search_pixels_x search_pixels_y] = meshgrid...
    (x-search_factor*cell_radius_x:x+search_factor*cell_radius_x, ...
    y-search_factor*cell_radius_y:y+search_factor*cell_radius_y);
% get rid of any pixels out of bounds
search_pixels_x = search_pixels_x(:);
search_pixels_y = search_pixels_y(:);
oob = false(size(search_pixels_x));
oob = search_pixels_x <= 0 | search_pixels_y <= 0 |...
    search_pixels_x > M | search_pixels_y > N;
search_pixels_x(oob) = [];
search_pixels_y(oob) = [];

% correlate by the third power of the trace, emphasizes peaks
search_pixels_indx = sub2ind([N M],search_pixels_y(:),search_pixels_x(:));
corrcoef_grid = zeros(N,M);
for i = search_pixels_indx';
    temp_corr = corrcoef(curr_trace.^3, ...
        double(handles.im(i,:)).^3);
    corrcoef_grid(i) = temp_corr(2);
end

corrcoef_grid = corrcoef_grid.^3;

corrcoef_threshold = mean(corrcoef_grid(search_pixels_indx)) + ...
    1*nanstd(corrcoef_grid(search_pixels_indx));

binaryImg = [];
binaryImg = corrcoef_grid >= corrcoef_threshold; % threshold image
if sum(binaryImg(:)) < 0.5*cell_area
    disp('No cell detected')
    return
end
binaryCell = [];
binaryCell = bwselect(binaryImg,x,y); % select cell that was clicked
if sum(binaryCell(:)) < 0.5*cell_area
    disp('No coherent cell detected')
    return
end

first_nonzero = find(binaryCell > 0,1);
[y_nonzero x_nonzero] = ind2sub([N M],first_nonzero);
roi_edge = bwtraceboundary(binaryCell,[y_nonzero x_nonzero],'N'); % find boundary, in order

% Get difference in distance to center across border
% function from file exchange
num_verticies = 10;
roi_edge_downsample = reduce_poly(roi_edge',num_verticies)';

% Repeat first point last
roi_edge_downsample(end+1,:) = roi_edge_downsample(1,:);

% create ROI from detected cell
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1) + 1;

polygon.ROI{n_polygon} = fliplr(roi_edge_downsample); % make x,y coordinates
polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});

set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);

% update handles
guidata(hObject, handles);


% --- Executes on button press in adaptHelpROI.
function adaptHelpROI_Callback(hObject, eventdata, handles)
% hObject    handle to adaptHelpROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

roi_numbers = get(handles.roi_numbers,'value');
polygon = get(handles.roiList,'UserData');
axes(handles.guiAxes);

% get centers of polygons
roi_centers = zeros(length(polygon.ROI),2);
for i = 1:length(polygon.ROI);
    % get center x and y from polygon boundaries
    center_x = mean(polygon.ROI{i}(:,1));
    center_y = mean(polygon.ROI{i}(:,2));
    roi_centers(i,1) = center_x;
    roi_centers(i,2) = center_y;
    delete(polygon.handles{i});
end

% clear loaded ROIs
set(handles.roiList,'String','');
set(handles.roiList,'UserData',[]);
axes(handles.guiAxes);

% use those centers as new x,y values for helping
if ~isfield(handles,'im_smooth')
    disp 'NOT SMOOTHING AT THE MOMENT'
    M = size(handles.im,1);
    N = size(handles.im,2);
    P = size(handles.im,3);
    im_long = reshape(handles.im,M*N,P);
    im_smooth = im_long;
    %     im_smooth = zeros(size(handles.im));
    %     M = size(handles.im,1);
    %     N = size(handles.im,2);
    %     P = size(handles.im,3);
    %     im_long = reshape(handles.im,M*N,P);
    %     im_smooth = zeros(M*N,P);
    %     disp('Smoothing movie...')
    %     parfor i = 1:M*N;
    %         im_smooth(i,:) = smooth(im_long(i,:),10);
    %     end
    %     disp('done')
    handles.im_smooth = im_smooth;
    % update handles
    guidata(hObject, handles);
end

% loop through all ROIs, use centers as seed points for new ROIs
axes(handles.guiAxes);
disp('Adapting ROIs')
for cell_num = 1:size(roi_centers,1);
    x = round(roi_centers(cell_num,1));
    y = round(roi_centers(cell_num,2));
    ind = sub2ind(size(mean(handles.im,3)),y,x);
    corrcoef_grid = [];
    corrcoef_grid = [handles.im_smooth(ind,:)*handles.im_smooth']' ./ sum(handles.im_smooth,2);
    corrcoef_grid_reshape = reshape(corrcoef_grid,size(mean(handles.im,3)));
    % threshold the corrcoef_grid_reshape
    corrcoef_threshold = nanmean(corrcoef_grid) + 2*nanstd(corrcoef_grid);
    binaryImg = [];
    binaryImg = corrcoef_grid_reshape >= corrcoef_threshold; % threshold image
    if any(any(binaryImg))
        binaryCell = [];
        binaryCell = bwselect(binaryImg,x,y); % select cell that was clicked
        if any(any(binaryCell))
            first_nonzero = find(binaryCell > 0,1);
            [y_nonzero x_nonzero] = ind2sub(size(binaryCell),first_nonzero);
            roi_edge = bwtraceboundary(binaryCell,[y_nonzero x_nonzero],'N'); % find boundary, in order
            % create ROI from detected cell
            polygon = get(handles.roiList,'UserData');
            n_polygon = size(get(handles.roiList,'string'),1) + 1;
            
            polygon.ROI{n_polygon} = fliplr(roi_edge); % make x,y coordinates
            polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
            set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
            
            set(handles.roiList,'string',num2str((1:n_polygon)'));
            set(handles.roiList,'UserData',polygon);
            
            % update handles
            guidata(hObject, handles);
        end
    end
    disp(['Adapted ' num2str(cell_num) '/' num2str(size(roi_centers,1))]);
end
% update handles
guidata(hObject, handles);


% --- Executes on button press in maxImg.
function maxImg_Callback(hObject, eventdata, handles)
% hObject    handle to maxImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxImg

maxImg = get(handles.maxImg,'value');
tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

if maxImg == 0;
    set(handles.imgSlider,'Enable','On');
    tiff_fullfilename = get(handles.tiffFile,'string');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
    
elseif maxImg == 1;
    set(handles.imgSlider,'Enable','Off');
    axes(handles.guiAxes);
    im_h = findobj('type','image');
    delete(im_h(1));
    hold on;
    imagesc(reshape(max(handles.im,[],2),N,M),'CDataMapping','scaled');
    colormap(gray)
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on button press in showBehavior.
function showBehavior_Callback(hObject, eventdata, handles)
% hObject    handle to showBehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showBehavior

% this assumes that you're looking at the decimated entire movie, and that
% doesn't make corrections for unimaged frames. basically, just a way to
% get an idea of what the lever was doing at the general time.
showBehavior = get(handles.showBehavior,'value');

if showBehavior == 1
    
    animal  = input('Animal: ','s');
    day = input('Day: ','s');
    
    if ispc
        data_path = ['Z:\Data\ImagingRig3\' day filesep animal];
    else
        data_path = ['/usr/local/lab/Data/ImagingRig3/' day filesep animal];
    end
    
    % Get behavior
    % [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    bhv_file = strncmp('data_@',dir_filenames,6);
    bhv_filename = [cell2mat(dir_filenames(bhv_file == 1))];
    
    warning off
    bhv = load([data_path filesep bhv_filename],'-MAT');
    warning on
    parsed_events = bhv.saved_history.ProtocolsSection_parsed_events;
    
    num_trials = length(parsed_events);
    
    all_reward = cellfun(@(x) x.states.reward, parsed_events,'UniformOutput',0);
    all_cue = cellfun(@(x) x.states.cue, parsed_events,'UniformOutput',0);
    all_iti = cellfun(@(x) x.states.iti, parsed_events,'UniformOutput',0);
    all_licks = cellfun(@(x) x.pokes.C(:,1), parsed_events,'UniformOutput',0);
    
    % get trial starts and ends
    all_trial_start = cellfun(@(x) x.states.state_0(1,2), parsed_events,'UniformOutput',0);
    all_trial_end = cellfun(@(x) x.states.state_0(2,1), parsed_events,'UniformOutput',0);
    
    reward_trials = [cellfun(@(x) ~isempty(x), all_reward)];
    
    % NOTE ON SETUP: C = lickport, L = lever
    % Index if/which lever press ended the cue state
    lever_cue_indx = cellfun(@(x) find(x.pokes.L(:,1) == x.states.cue(:,2)), parsed_events,'UniformOutput',0);
    
    % Get lever force
    
    % find and load the XSG files
    acq_path = [data_path filesep 'AP' '00' animal(3:4)];
    dir_currfolder_acq = dir(acq_path);
    dir_acq_filenames = {dir_currfolder_acq.name};
    acq_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.xsg')), dir_acq_filenames);
    
    
    acq_filename_all = dir_acq_filenames(acq_filename_indx);
    % make sure they're in order
    acq_filename_all = sort(acq_filename_all);
    
    lever_force_split = cell(length(acq_filename_all),1);
    
    lever_force = [];
    trial_list = [];
    lever_force = [];
    bitcode_trace = [];
    trial_stitch = [];
    % At the moment, concatenate all and pretend no interruption
    disp('ASSUMING: 8 FULL FILES FOR EACH LOOP')
    
    % Get framerate
    
    % Woah woah WOAH woah woah - this isn't always correct! This depends on
    % hitting the 'measure framerate' button at the beginning, and even then
    % it may not be right. Here's what we'll do instead: above, get the
    % framerate by dividing the number of frames by the length of the xsg in
    % seconds. Note that this is not robust to dropped frames.
    % Edit: at the moment, I know how many frames there should be, so I'll just
    % hardcode it.
    % 2nd Edit: This doesn't make sense - the XSG is slightly longer than the
    % acq, the ONLY information is 'framerate' in the header, which is really
    % bullshit that it's not. For now, assume it's right, but fix for AP60
    % which I know didn't have measured framerate for the first couple days
    % disp('ASSUMING 4000 FRAMES per XSG!');
    % framerate = 4000/(mean(cellfun(@length,lever_force_split))/xsg_sample_rate);
    
    % [img_filename img_path] = uigetfile('*.tif','Pick sample file','Multiselect','off');
    dir_currfolder = dir(data_path);
    dir_filenames = {dir_currfolder.name};
    tiff_filename_indx = cellfun(@(x) ~isempty(strfind(x,'.tif')), dir_filenames);
    
    img_filename = dir_filenames(tiff_filename_indx);
    
    img_info = imfinfo([data_path filesep img_filename{1}]);
    img_info = img_info(1).ImageDescription;
    [img_parameter img_value] = strread(img_info,'%s %s', 'delimiter','=\n');
    
    framerate_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.scanFrameRate')),img_parameter,'UniformOutput',0));
    framerate = str2num(img_value{framerate_indx});
    
    if strcmp(img_filename{1}(8:11),'AP60')
        framerate = 28.7224;
        disp('OVERRIDING FRAMERATE FOR AP60! NOT PERMANENT SOLUTION TO THIS PROBLEM!')
    end
    
    for i = 1:length(acq_filename_all);
        xsg = load([acq_path filesep acq_filename_all{i}],'-MAT');
        
        % Get all trial offsets in samples
        xsg_sample_rate = xsg.header.acquirer.acquirer.sampleRate;
        
        % Get trials in raw samples since started
        curr_trial_list = [];
        curr_trial_list = AP_ReadBitCode([acq_path filesep acq_filename_all{i}]);
        if ~isempty(curr_trial_list)
            curr_trial_list(:,1) = curr_trial_list(:,1) + length(lever_force)/xsg_sample_rate;
            trial_list = [trial_list; curr_trial_list];
        end
        
        % ignore this part for just general display purposes
        
        % Cut out pieces of the lever trace that weren't imaged (dropped
        % frames) - Assume they were clipped at the end
        %         curr_loop_frames = sum(cellfun(@(x) size(x,2),roi_trace_long_split((8*(i-1)+1):8*i)));
        %         curr_imaged_lever_samples = round((curr_loop_frames/framerate)*xsg_sample_rate);
        %         curr_lever_trace = xsg.data.acquirer.trace_2(1:curr_imaged_lever_samples);
        curr_lever_trace = xsg.data.acquirer.trace_2;
        
        bitcode_trace = [bitcode_trace; xsg.data.acquirer.trace_1];
        lever_force_split{i} = xsg.data.acquirer.trace_2;
        lever_force = [lever_force; curr_lever_trace];
        trial_stitch = [trial_stitch trial_list(end,2)];
        disp(['Loaded xsg ' num2str(i) '/' num2str(length(acq_filename_all))]);
    end
    
    % don't resample, just show as is
    
    %     % resample lever force to be 1:1 with numframes
    %     numframes = size(roi_trace_long,2);
    %     % to make resampling work: numframes must be even number
    %     % if it's not, cut out one frame at the end
    %     if mod(numframes,2) ~= 0
    %         roi_trace_long(:,end) = [];
    %         numframes = size(roi_trace_long,2);
    %     end
    % resample lever force to so # samples = # frames
    %     [n d] = rat(numframes/length(lever_force));
    %     lever_force_resample = resample(lever_force,n,d);
    % add or delete last sample if necessary
    %     if length(lever_force_resample) < numframes
    %         lever_force_resample(end+1:numframes) = 0;
    %     elseif length(lever_force_resample) > numframes
    %         lever_force_resample = lever_force_resample(1:numframes);
    %     end
    
    % resample it so that one sample ~ one frame
    [n d] = rat(framerate/xsg_sample_rate);
    lever_force_resample = resample(lever_force,n,d);
    
    % Get trial offsets (seconds)
    xsg_bhv_offset = nan(num_trials,1);
    for curr_trial = 1:num_trials
        xsg_trial_indx = find(trial_list(:,2) == curr_trial,1);
        if ~isempty(xsg_trial_indx)
            xsg_bhv_offset(curr_trial) = ...
                trial_list(xsg_trial_indx,1) - ...
                all_trial_start{curr_trial}(1);
        else
            continue;
        end
    end
    
    % Get dispatcher times in frames
    cue_frames = [];
    reward_frames = [];
    lick_frames = [];
    for curr_trial = 1:length(xsg_bhv_offset);
        if isnan(xsg_bhv_offset(curr_trial))
            continue
        end
        
        % relative cue time
        cue_frames = [cue_frames;(all_cue{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        % relative reward time
        if ~isempty(all_reward{curr_trial})
            reward_frames = [reward_frames;(all_reward{curr_trial}(1,1) + xsg_bhv_offset(curr_trial))*framerate];
        end
        % relative lick time
        lick_frames = [lick_frames;(all_licks{curr_trial} + xsg_bhv_offset(curr_trial))*framerate];
    end
    
    
    axes(handles.behaviorPlot)
    hold on
    plot(lever_force_resample,'k')
    for i = 1:length(reward_frames)
        line([reward_frames(i) reward_frames(i)],ylim,'linestyle','--',...
            'color','g');
    end
    for i = 1:length(cue_frames)
        line([cue_frames(i) cue_frames(i)],ylim,'linestyle','--',...
            'color','b');
    end
    
    % save the resampled lever force
    handles.lever_force_resample = lever_force_resample;
    % update handles
    guidata(hObject, handles);
elseif showBehavior == 0;
    axes(handles.behaviorPlot)
    clear axes
end



function cell_diameter_Callback(hObject, eventdata, handles)
% hObject    handle to cell_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_diameter as text
%        str2double(get(hObject,'String')) returns contents of cell_diameter as a double


% --- Executes during object creation, after setting all properties.
function cell_diameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_diameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lum_radius_Callback(hObject, eventdata, handles)
% hObject    handle to lum_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lum_radius as text
%        str2double(get(hObject,'String')) returns contents of lum_radius as a double


% --- Executes during object creation, after setting all properties.
function lum_radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lum_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function blur_std_Callback(hObject, eventdata, handles)
% hObject    handle to blur_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blur_std as text
%        str2double(get(hObject,'String')) returns contents of blur_std as a double


% --- Executes during object creation, after setting all properties.
function blur_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blur_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lumnorm.
function lumnorm_Callback(hObject, eventdata, handles)
% hObject    handle to lumnorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Normalizing luminance...')

lum_radius = str2num(get(handles.lum_radius,'string'));

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
numframes = size(handles.im,2);
% reshape to regular dimension movie
% warning: committing yourself to double memory here
handles.im = single(reshape(handles.im,N,M,numframes));

% make the luminance normalization filter, filter the avg
h = fspecial('average', lum_radius);
im_avg = mean(handles.im,3);
im_avg_smooth = imfilter(im_avg,h);
% smooth/lumnorm image with gaussian filter
disp('Finished (%):')
fprintf('%2d',0)
for i = 1:numframes
    %handles.im(:,:,i) = imfilter((handles.im(:,:,i)-im_avg)./(im_avg+im_avg_smooth),h2).^3;%(handles.im(:,:,i)-im_avg)./(im_avg*im_avg_smooth)
    handles.im(:,:,i) = (handles.im(:,:,i)-im_avg)./(im_avg+im_avg_smooth);
    fprintf('%c%c%2d',8,8,round(100*i/numframes));
end
% put back to being frames in rows
handles.im = reshape(handles.im,N*M,numframes);
% update figure
axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(reshape(handles.im(:,1),N,M),'CDataMapping','scaled');
xlim([0 M]);
ylim([0 N]);
colormap(gray)
imgMin = min(handles.im(:));
imgMax = max(handles.im(:));
set(handles.imgMin,'string',0);
set(handles.imgMax,'string',imgMax/3);
caxis([0 imgMax/3]);
numframes=size(handles.im,2);
set(handles.imgSlider,'Min',1);
if numframes ~= 1;
    set(handles.imgSlider,'Enable','on');
    set(handles.imgSlider,'Max',numframes);
    set(handles.imgSlider,'Value',1);
    set(handles.imgSlider,'SliderStep',[1/(numframes-1), 10/(numframes-1)]);
else
    set(handles.imgSlider,'Enable','off');
end

% flip the ROIs to be on top
set(gca,'Children',flipud(get(gca,'Children')));

% update handles
guidata(hObject, handles);

disp('Done.')


% --- Executes on button press in load_imagej.
function load_imagej_Callback(hObject, eventdata, handles)
% hObject    handle to load_imagej (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
polygon = get(handles.roiList,'UserData');

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);

[roi_file roi_path] = uigetfile('.zip','Pick ImageJ ROI zip',tiff_filename);

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

axes(handles.guiAxes);

% % delete existing polygons and lines
% currpoly = findobj(handles.guiAxes,'tag','impoly');
% delete(currpoly)
% currpoly = findobj(handles.guiAxes,'type','line');
% delete(currpoly)

imagej_polygon = readimagejroi([roi_path roi_file]);
curr_polygons = length(polygon.ROI);
for n_polygon = 1:length(imagej_polygon)
    polygon.ROI{curr_polygons+n_polygon} = imagej_polygon{n_polygon}.mnCoordinates;
    polygon.ROI{curr_polygons+n_polygon} = [polygon.ROI{curr_polygons+n_polygon};polygon.ROI{curr_polygons+n_polygon}(1,:)];
    polygon.handles{curr_polygons+n_polygon} = line(polygon.ROI{curr_polygons+n_polygon}(:,1),polygon.ROI{curr_polygons+n_polygon}(:,2));
    set(polygon.handles{curr_polygons+n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
end

set(handles.roiList,'string',num2str((1:curr_polygons+n_polygon)'));
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in load_zebraroi.
function load_zebraroi_Callback(hObject, eventdata, handles)
% hObject    handle to load_zebraroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tiff_fullfilename = get(handles.tiffFile,'string');
[zebraroi_file zebraroi_path] = uigetfile('.zebraroi','Pick AutoZebra file');
load([zebraroi_path zebraroi_file],'-MAT');

imageinfo=imfinfo(tiff_fullfilename,'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

axes(handles.guiAxes);

% get rid of any empty polygons
zebra_roi_indx = cellfun(@(x) ~isempty(x),zebra_roi);
zebra_roi = zebra_roi(zebra_roi_indx);

% create ROI outlines from filled in ROIs
for n_polygon = 1:length(zebra_roi)
    binaryCell = false(N,M);
    binaryCell(sub2ind([N M],zebra_roi{n_polygon}(:,1),zebra_roi{n_polygon}(:,2))) = true;
    
    %first_nonzero = find(binaryCell > 0,1);
    %[y_nonzero x_nonzero] = ind2sub([N M],first_nonzero);
    roi_edge = bwtraceboundary(binaryCell,[zebra_roi{n_polygon}(1,1) zebra_roi{n_polygon}(1,2)],'N'); % find boundary, in order
    % create ROI from detected cell
    polygon = get(handles.roiList,'UserData');
    n_polygon = size(get(handles.roiList,'string'),1) + 1;
    
    polygon.ROI{n_polygon} = roi_edge;%fliplr(roi_edge); % make x,y coordinates
    polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1), polygon.ROI{n_polygon}(:,2));
    set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
    
    set(handles.roiList,'string',num2str((1:n_polygon)'));
    set(handles.roiList,'UserData',polygon);
    
    % update handles
    guidata(hObject, handles);
end

set(handles.roiList,'string',num2str((1:n_polygon)'));
set(handles.roiList,'UserData',polygon);


% --- Executes on button press in show_zebra_roi.
function show_zebra_roi_Callback(hObject, eventdata, handles)
% hObject    handle to show_zebra_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Normalizing luminance...')

lum_radius = str2num(get(handles.lum_radius,'string'));

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
numframes = size(handles.im,2);
% reshape to regular dimension movie
handles.im = reshape(handles.im,N,M,numframes);

% make the luminance normalization filter, filter the avg
h = fspecial('disk', lum_radius);
im_avg = mean(handles.im,3);
im_avg_smooth = double(imfilter(im_avg,h));
% smooth/lumnorm image with gaussian filter
disp('Finished (%):')
fprintf('%2d',0)
im_zebraroi = zeros(size(im_avg),'double');
for i = 1:numframes
    curr_frame = zeros(size(im_avg),'double');
    curr_frame = (double(handles.im(:,:,i))-double(im_avg))./(double(im_avg)+double(im_avg_smooth));
    curr_frame = curr_frame.^3;
    curr_frame = imfilter(curr_frame,h);
    im_zebraroi = im_zebraroi+curr_frame./numframes;
    fprintf('%c%c%2d',8,8,round(100*i/numframes));
end

% update figure
% axes(handles.guiAxes);
% set(handles.imgSlider,'Enable','Off');
% axes(handles.guiAxes);
% im_h = findobj('type','image');
% delete(im_h(1));
figure;
imagesc(im_zebraroi);
xlim([0 M]);
ylim([0 N]);
colormap(gray)
imgMin = min(im_zebraroi(:));
imgMax = max(im_zebraroi(:));
caxis([0 imgMax/30]);

% update handles
handles.im = reshape(handles.im,N*M,numframes);

guidata(hObject, handles);

disp('Done.')


% --- Executes on button press in roi_label.
function roi_label_Callback(hObject, eventdata, handles)
% hObject    handle to roi_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in roi_labelmaker.
function roi_labelmaker_Callback(hObject, eventdata, handles)
% hObject    handle to roi_labelmaker (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.roi_labelmaker = AP_roi_labelmaker('AP_cellROI_pro',hObject);
guidata(hObject,handles);



function bg_roi_outer_Callback(hObject, eventdata, handles)
% hObject    handle to bg_roi_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bg_roi_outer as text
%        str2double(get(hObject,'String')) returns contents of bg_roi_outer as a double


% --- Executes during object creation, after setting all properties.
function bg_roi_outer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bg_roi_outer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bg_roi_inner_Callback(hObject, eventdata, handles)
% hObject    handle to bg_roi_inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bg_roi_inner as text
%        str2double(get(hObject,'String')) returns contents of bg_roi_inner as a double


% --- Executes during object creation, after setting all properties.
function bg_roi_inner_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bg_roi_inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bg_rois.
function bg_rois_Callback(hObject, eventdata, handles)
% hObject    handle to bg_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bg_rois

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
aspect_ratio = M/N;

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1);
bg_roi_inner = str2num(get(handles.bg_roi_inner,'String'));
bg_roi_outer = str2num(get(handles.bg_roi_outer,'String'));
bg_roi_outer_full = bg_roi_outer + bg_roi_inner;

% if there are already drawn background ROIs, delete their lines
if isfield(polygon,'bghandles')
    delete(polygon.bghandles{:});
end
    
% loop through ROIs, draw background ROIs with given inner/outer radius
for curr_roi = 1:n_polygon
    
    % get center of current polygon
    [area roi_center_x roi_center_y] = ...
        polycenter(polygon.ROI{curr_roi}(:,1), ...
        polygon.ROI{curr_roi}(:,2));
    
    % center current polygon at zero
    temp_roi = [polygon.ROI{curr_roi}(:,1) - roi_center_x ...
        polygon.ROI{curr_roi}(:,2) - roi_center_y];
    
%     %%%
%     % TEMP: MAKE INNER/OUTER'
%     %%%
%     
%     % expand polygon by inner and outer diameters
%     temp_roi_inner = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_inner ...
%         temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_inner/aspect_ratio];   
%     temp_roi_outer = [temp_roi(:,1) + sign(temp_roi(:,1))*round(bg_roi_outer_full/2) ...
%         temp_roi(:,2) + sign(temp_roi(:,2))*round(bg_roi_outer_full/2)/aspect_ratio];
%     
%     % connect the inner and outer sections
%     bg_roi = [temp_roi_outer;temp_roi_inner;temp_roi_outer(1,:)];
%     bg_roi_centered_1 = [bg_roi(:,1) + roi_center_x ...
%         bg_roi(:,2) + roi_center_y];
%     
%         % expand polygon by inner and outer diameters
%     temp_roi_inner = [temp_roi(:,1) + sign(temp_roi(:,1))*round(bg_roi_outer_full/2) ...
%         temp_roi(:,2) + sign(temp_roi(:,2))*round(bg_roi_outer_full/2)/aspect_ratio];   
%     temp_roi_outer = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_outer_full ...
%         temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_outer_full/aspect_ratio];
%     
%     % connect the inner and outer sections
%     bg_roi = [temp_roi_outer;temp_roi_inner;temp_roi_outer(1,:)];
%     bg_roi_centered_2 = [bg_roi(:,1) + roi_center_x ...
%         bg_roi(:,2) + roi_center_y];
    
        % expand polygon by inner and outer diameters
    temp_roi_inner = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_inner ...
        temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_inner/aspect_ratio];   
    temp_roi_outer = [temp_roi(:,1) + sign(temp_roi(:,1))*bg_roi_outer_full ...
        temp_roi(:,2) + sign(temp_roi(:,2))*bg_roi_outer_full/aspect_ratio];
    
    % connect the inner and outer sections
    bg_roi = [temp_roi_outer;temp_roi_inner;temp_roi_outer(1,:)];
    bg_roi_centered = [bg_roi(:,1) + roi_center_x ...
        bg_roi(:,2) + roi_center_y];
    
    %%%
    %
    %%%
    
    % save the bgROI
    polygon.bgROI{curr_roi} = bg_roi_centered;
    
    % eliminate uncorrelated pixels within background ROI to avoid
    % subtracting non-contaminating information
    
    % get ROI trace
    curr_roi_mask = poly2mask(polygon.ROI{curr_roi}(:,1)',...
        polygon.ROI{curr_roi}(:,2)',N,M);
    curr_roi_trace = mean(handles.im(curr_roi_mask(:),:));
    
    % get background trace
    curr_bgroi_mask = poly2mask(polygon.bgROI{curr_roi}(:,1)',...
        polygon.bgROI{curr_roi}(:,2)',N,M);
    curr_bgroi_pixel_trace = handles.im(curr_bgroi_mask(:),:);
    
%     % split this in random halves for ICA ::: this gave same result?
%     num_bgroi_px = sum(curr_bgroi_mask(:));
%     bgroi_mask_indx = find(curr_bgroi_mask);
%     bgroi_rand_indx = randperm(num_bgroi_px);
%     bgroi_mask_half_1 = bgroi_mask_indx(bgroi_rand_indx(1:floor(num_bgroi_px/2)));
%     bgroi_mask_half_2 = bgroi_mask_indx(bgroi_rand_indx(floor(num_bgroi_px/2)+1:end));
%     
%     bgroi_trace_1 = mean(handles.im(bgroi_mask_half_1,:));
%     bgroi_trace_2 = mean(handles.im(bgroi_mask_half_2,:));
    
%     % get correlation between background ROI and drawn ROI
%     % correlate by the third power of the trace, emphasizes peaks
%     corrcoef_grid_roi = zeros(N,M);
%     for i = find(curr_roi_mask)';
%         temp_corr = corrcoef(curr_roi_trace, ...
%             double(handles.im(i,:)));
%         corrcoef_grid_roi(i) = temp_corr(2);
%     end
%     
%     corrcoef_grid_bgroi = zeros(N,M);
%     for i = find(curr_bgroi_mask)';
%         temp_corr = corrcoef(curr_roi_trace, ...
%             double(handles.im(i,:)));
%         corrcoef_grid_bgroi(i) = temp_corr(2);
%     end
    
    % draw the background ROIs
    polygon.bghandles{curr_roi} = line(polygon.bgROI{curr_roi}(:,1), ...
        polygon.bgROI{curr_roi}(:,2),'color','m','linestyle','--');
    
    % save the masks for this too: easier to manupilate later 
    polygon.bgROI_mask{curr_roi} = curr_bgroi_mask;
end

set(handles.roiList,'UserData',polygon);


% --- Executes on button press in fix_bg_rois.
function fix_bg_rois_Callback(hObject, eventdata, handles)
% hObject    handle to fix_bg_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bg_rois_Callback(hObject, eventdata, handles);

tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;
aspect_ratio = M/N;

axes(handles.guiAxes);
polygon = get(handles.roiList,'UserData');
n_polygon = size(get(handles.roiList,'string'),1);
bg_cutoff = str2num(get(handles.bg_cutoff,'string'));

% loop through bgROIs, find pixels which contain transients that exceed a
% standard difference from the main ROI
for curr_roi = 1:n_polygon
     
    % get ROI trace
    curr_roi_mask = poly2mask(polygon.ROI{curr_roi}(:,1)',...
        polygon.ROI{curr_roi}(:,2)',N,M);
    curr_roi_trace = mean(handles.im(curr_roi_mask(:),:));
    
    % get background trace
    curr_bgroi_pixel_trace = handles.im(polygon.bgROI_mask{curr_roi}(:),:);
    
    bg_diff = ...
        bsxfun(@minus,double(curr_bgroi_pixel_trace),double(curr_roi_trace));
    
    % define outliers by the distribution of differences. There is a
    % gaussian (sometimes a little buried) - so just estimate the mean with
    % mode and the std
    
    % get mode rounded to nearest integer
    bg_diff_mode = mode(round(bg_diff(:)*1))/1;
    bg_diff_std = std(bg_diff(:));
    bg_cutoff = bg_diff_mode + bg_cutoff*bg_diff_std;
    
    bg_trace_outliers = ...
        bsxfun(@gt,curr_bgroi_pixel_trace,curr_roi_trace+bg_cutoff);
    bg_mask_indx = find(polygon.bgROI_mask{curr_roi});
    bg_pixel_outliers = bg_mask_indx(any(bg_trace_outliers,2));
    
    polygon.bgROI_mask{curr_roi}(bg_pixel_outliers) = false;
    
    % if there are < 10 pixels, use whole bgroi
    if sum(polygon.bgROI_mask{curr_roi}(:)) < 10
        polygon.bgROI_mask{curr_roi}(bg_pixel_outliers) = true;
        %disp(['Warning: bgROI has less than 10 pixels, using whole bgROI: ' num2str(curr_roi)]);
        set(polygon.bghandles{curr_roi},'color','r');
        
    % give a warning if < 50 pixels
    elseif sum(polygon.bgROI_mask{curr_roi}(:)) < 50
        %disp(['Warning: bgROI has less than 50 pixels: ' num2str(curr_roi)]);
        set(polygon.bghandles{curr_roi},'color',[1 0.5 0]);
        
    else
        set(polygon.bghandles{curr_roi},'color','g');  
        
    end
    
    
    
end

set(handles.roiList,'UserData',polygon);


% --- Executes on button press in fine_correct.
function fine_correct_Callback(hObject, eventdata, handles)
% hObject    handle to fine_correct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fine_correct


% --- Executes on button press in show_bgrois.
function show_bgrois_Callback(hObject, eventdata, handles)
% hObject    handle to show_bgrois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_bgrois

show_bgrois = get(handles.show_bgrois,'value');
polygon = get(handles.roiList,'UserData');

bgroi_display = max(cat(3,polygon.bgROI_mask{:}),[],3);



function bg_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to bg_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bg_cutoff as text
%        str2double(get(hObject,'String')) returns contents of bg_cutoff as a double


% --- Executes during object creation, after setting all properties.
function bg_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bg_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roi_auto_align_batch.
function roi_auto_align_batch_Callback(hObject, eventdata, handles)
% hObject    handle to roi_auto_align_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is to align ROIs to multiple selected summed movies. Select the
% template ROI and all tiff files (the template ROI must come from one of
% them and contain the full title of the movie it is derived from)

% Pick template ROI file and select movie
[roi_file, roi_path] = uigetfile('*.roi','Choose template ROI');
[tiff_files, tiff_path] = uigetfile('*.tif','Choose tiff files','multiselect','on');
tiff_files = sort(tiff_files);

% Determine which movie the current ROI was derived from
roi_tiff = cellfun(@(x) any(strfind(roi_file,x(1:end-4))),tiff_files);

% Set up saved variables
all_avg_ims = cell(size(tiff_files));
all_max_ims = cell(size(tiff_files));
all_aligned_polygon = cell(size(tiff_files));

% Get the average image from the template file
waitbar_h = waitbar(0,'Creating average from reference day');

orig_file = tiff_files{roi_tiff};
imageinfo_orig = imfinfo([tiff_path orig_file]);
N_orig = imageinfo_orig(1).Height;
M_orig = imageinfo_orig(1).Width;
im_orig_avg = zeros(N_orig,M_orig);
im_orig_max = zeros(N_orig,M_orig);
for curr_frame = 1:length(imageinfo_orig)
    im_orig = double(imread([tiff_path orig_file],'tiff',curr_frame));
    im_orig_avg = im_orig_avg + im_orig/length(imageinfo_orig);
    im_orig_max = max(cat(3,im_orig_max,im_orig),[],3);
end
all_avg_ims{roi_tiff} = im_orig_avg;
all_max_ims{roi_tiff} = im_orig_max;

% Load in the template polygon 
template = load([roi_path roi_file],'-mat');
all_aligned_polygon{roi_tiff} = template.polygon.ROI;

% Loop through each file, auto-align ROIs (align each one to the template
% in order to not propogate errors)
waitbar(0,waitbar_h,'Auto-aligning template ROI')
for curr_file = find(~roi_tiff);

    % Don't bother trying to fine correct at the moment, removed that here
    
    [tiff_fullfilename] = [tiff_path tiff_files{curr_file}];
    imageinfo = imfinfo(tiff_fullfilename,'tiff');
    numframes=length(imageinfo);
    M = imageinfo(1).Width;
    N = imageinfo(1).Height;
    
    aspect_ratio = M/N;
    
    im_curr_avg = zeros(N,M);
    im_curr_max = zeros(N,M);
    for loadframe = 1:numframes
        curr_frame = double(imread(tiff_fullfilename,'tiff',loadframe,'Info',imageinfo));
        im_curr_avg = im_curr_avg + curr_frame./numframes;
        im_curr_max = max(cat(3,im_curr_max,curr_frame),[],3);
    end
    all_avg_ims{curr_file} = im_curr_avg;
    all_max_ims{curr_file} = im_curr_max;
        
    % Estimate average ROI radius
    polygon_areas = zeros(length(template.polygon.ROI),1);
    for curr_roi = 1:length(template.polygon.ROI)
        polygon_areas(curr_roi) = polyarea( ...
            template.polygon.ROI{curr_roi}(:,1),aspect_ratio*template.polygon.ROI{curr_roi}(:,2));
    end
    
    cell_radius = round(sqrt(mean(polygon_areas)/pi));
    
    % Luminance normalize average images
    h = fspecial('average', [cell_radius ceil(cell_radius/aspect_ratio)]);
    im_orig_avg_lumnorm = im_orig_avg./imfilter(im_orig_avg,h);
    im_curr_avg_lumnorm = im_curr_avg./imfilter(im_curr_avg,h);
    
    % Occasionally get NaN/Inf around edges, set to zero
    im_orig_avg_lumnorm(isinf(im_orig_avg_lumnorm) | ...
        isnan(im_orig_avg_lumnorm)) = 0;
    im_curr_avg_lumnorm(isinf(im_curr_avg_lumnorm) | ...
        isnan(im_curr_avg_lumnorm)) = 0;
    
    % Go through each ROI, figure out the best new center from the old center
    % this assumes little rotation and no movement of rois from straight load
    s = 5;
    s_aspect = round(s/aspect_ratio);
    
    aligned_polygon = template.polygon.ROI;
    uncorrected_flag = false(length(aligned_polygon),1);
    corr_offset_all = zeros(length(template.polygon.ROI),2);
    
    % estimate centers of all ROIs
    roi_centers = zeros(length(template.polygon.ROI),2);
    for curr_roi = 1:length(template.polygon.ROI)
        % estimate center of roi
        center_x = round(mean(template.polygon.ROI{curr_roi}(:,1)));
        center_y = round(mean(template.polygon.ROI{curr_roi}(:,2)));
        roi_centers(curr_roi,:) = [center_x center_y];
    end
    
    for curr_roi = 1:length(template.polygon.ROI)
        % estimate center of roi
        center_x = roi_centers(curr_roi,1);
        center_y = roi_centers(curr_roi,2);
        roi_snip = [];
        if center_y-s_aspect*cell_radius <= 0 || center_x-s*cell_radius <= 0 || ...
                center_y+s_aspect*cell_radius > N || center_x+s*cell_radius > M
            uncorrected_flag(curr_roi) = true;
            continue
        end
        % pull out s*cell_radius around the cell, cross correlate
        roi_snip = im_orig_avg_lumnorm(center_y-s_aspect*cell_radius:center_y+s_aspect*cell_radius, ...
            center_x-s*cell_radius:center_x+s*cell_radius);
        cc = normxcorr2(roi_snip,im_curr_avg_lumnorm);
        
        % get the point of max correlation
        xcorr_midpoint = round(length(cc)/2);
        [max_c,imax] = max(cc(:));
        [ypeak xpeak] = ind2sub(size(cc),imax(1));
        corr_offset = [(xpeak-s*cell_radius)-center_x ...
            (ypeak-s_aspect*cell_radius)-center_y];
        
        corr_offset_all(curr_roi,:) = corr_offset;
        
        % record the new positions based on center offset
        aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
            repmat(corr_offset,size(aligned_polygon{curr_roi},1),1);
    end
    
    % find abnormal offsets, reset them and flag them as uncorrected
    corr_offset_all_dist = sqrt(corr_offset_all(:,1).^2 + ...
        (s_aspect*corr_offset_all(:,2)).^2);
    
    % use 4x the distance of spread of middle 80% as offset cutoff
    offset_prctile = prctile(corr_offset_all_dist(corr_offset_all_dist ~= 0),[10 90]);
    offset_cutoff = 4*diff(offset_prctile);
    
    abnormal_offsets = find(corr_offset_all_dist < ...
        median(corr_offset_all_dist) - offset_cutoff | ...
        corr_offset_all_dist > ...
        median(corr_offset_all_dist) + offset_cutoff);
    aligned_polygon(abnormal_offsets) = template.polygon.ROI(abnormal_offsets);
    uncorrected_flag(abnormal_offsets) = true;
    
    % Estimate offset for unaligned ROIs
    median_corr_offset = median(corr_offset_all(~uncorrected_flag));
    for curr_roi = find(uncorrected_flag)'
        % estimate offset based on nearby ROIs, get distance
        curr_roi_center = [];
        curr_roi_center = repmat(roi_centers(curr_roi,:),length(roi_centers),1);
        curr_roi_dist = [];
        curr_roi_dist = sqrt(sum((roi_centers - curr_roi_center).^2,2));
        curr_roi_dist = [curr_roi_dist (1:length(curr_roi_dist))'];
        % don't use ROIs that weren't aligned
        curr_roi_dist = curr_roi_dist(~uncorrected_flag,:);
        % sort them by distance, use the top 3
        curr_roi_dist = sortrows(curr_roi_dist,1);
        
        corr_offset = median(corr_offset_all(curr_roi_dist(1:3,2),:));
        
        corr_offset_all(curr_roi,:) = corr_offset;
        
        aligned_polygon{curr_roi} = aligned_polygon{curr_roi} + ...
            repmat(corr_offset,size(aligned_polygon{curr_roi},1),1);
    end
  
    all_aligned_polygon{curr_file} = aligned_polygon;
    
    waitbar(curr_file/length(tiff_files),waitbar_h,'Auto-aligning template ROI')
end
close(waitbar_h);

% Plot all the aligned ROIs
sq_fig = ceil(sqrt(length(tiff_files)));
figure;
for curr_file = 1:length(tiff_files)
    subplot(sq_fig,sq_fig,curr_file);
    imagesc(all_avg_ims{curr_file});
    colormap(gray);caxis([0 2000]);

    % plot the old ROIs in dark gray
    for curr_roi = 1:length(template.polygon.ROI)
        roi1(curr_roi) = line(template.polygon.ROI{curr_roi}(:,1), ...
            template.polygon.ROI{curr_roi}(:,2),'color',[1 0 0]);
    end
    
    % plot the new ROIs in light gray
    for curr_roi = 1:length(all_aligned_polygon{curr_file})
        roi2(curr_roi) = line(all_aligned_polygon{curr_file}{curr_roi}(:,1), ...
            all_aligned_polygon{curr_file}{curr_roi}(:,2),'color',[0 1 0]);
    end
   
    % draw lines between old and new ROIs, label ROI numbers
    for curr_roi = 1:length(template.polygon.ROI)
        center_x1 = round(mean(template.polygon.ROI{curr_roi}(:,1)));
        center_y1 = round(mean(template.polygon.ROI{curr_roi}(:,2)));
        
        center_x2 = round(mean(all_aligned_polygon{curr_file}{curr_roi}(:,1)));
        center_y2 = round(mean(all_aligned_polygon{curr_file}{curr_roi}(:,2)));
        
        line([center_x1 center_x2],[center_y1 center_y2], ...
            'color','w','linewidth',2)
        
        roi_numbers_h(curr_roi) = text(center_x1,center_y1,num2str(curr_roi), ...
        'color','white','FontWeight','bold','FontSize',12);
    end
      
    title(tiff_files{curr_file},'Interpreter','none');
    
    % Save the aligned template ROI
    clear polygon
    polygon.ROI = all_aligned_polygon{curr_file};
    save([roi_path tiff_files{curr_file}(1:end-4) '_template.roi'],'polygon');
    
end


% --- Executes on button press in compare_rois.
function compare_rois_Callback(hObject, eventdata, handles)
% hObject    handle to compare_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is to manually compare all ROIs across all sessions, with the option
% to move the roi, change the verticies, or delete (can be opened with this
% GUI or without, they're independent)

compare_script = menu('Choose ROI compare format','Gray','Color','Gray - known events');
if compare_script == 1
    AP_cellROI_compareROI
elseif compare_script == 2
    AP_cellROI_compareROI_color;
elseif compare_script == 3;
    AP_cellROI_compareROI_known_events
end



% --- Executes on button press in dendrite_rois.
function dendrite_rois_Callback(hObject, eventdata, handles)
% hObject    handle to dendrite_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AP_cellROI_dendriteROI;



function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channels as text
%        str2double(get(hObject,'String')) returns contents of channels as a double


% --- Executes during object creation, after setting all properties.
function channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in compare_rois_field.
function compare_rois_field_Callback(hObject, eventdata, handles)
% hObject    handle to compare_rois_field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AP_compare_rois_field;



