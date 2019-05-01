function varargout = AP_cellROI(varargin)
% AP_CELLROI MATLAB code for AP_cellROI.fig
%      AP_CELLROI, by itself, creates a new AP_CELLROI or raises the existing
%      singleton*.
%
%      H = AP_CELLROI returns the handle to a new AP_CELLROI or the handle to
%      the existing singleton*.
%
%      AP_CELLROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AP_CELLROI.M with the given input arguments.
%
%      AP_CELLROI('Property','Value',...) creates a new AP_CELLROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AP_cellROI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AP_cellROI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AP_cellROI

% Last Modified by GUIDE v2.5 31-Jan-2012 16:37:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AP_cellROI_OpeningFcn, ...
    'gui_OutputFcn',  @AP_cellROI_OutputFcn, ...
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


% --- Executes just before AP_cellROI is made visible.
function AP_cellROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AP_cellROI (see VARARGIN)

% Choose default command line output for AP_cellROI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AP_cellROI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AP_cellROI_OutputFcn(hObject, eventdata, handles)
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
AP_cellROI_createROI

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

if isfield(handles,'autosortImg')
    handles = rmfield(handles,'autosortImg');
    handles = rmfield(handles,'guiAxes');
    guidata(gcbo,handles);
end

% use those centers as new x,y values for helping
if isfield(handles,'im_smooth')
    handles = rmfield(handles,'im_smooth');
    % update handles
    guidata(hObject, handles);
end

[tiff_filename,tiff_path]=uigetfile('*.tif','Choose File','Multiselect','off');
set(handles.tiffFile,'string',[tiff_path tiff_filename]);
cd(tiff_path);
if isfield(handles,'guiAxes') == 0;
    figure('KeyPressFcn', {@keyPress, hObject}) % set up figure for keyboard shortcuts
    handles.guiAxes = axes;
    guidata(gcbo,handles);
end

%%%% load in image file

imageinfo=imfinfo([tiff_path tiff_filename],'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;
clear imageinfo

im = zeros(N,M,numframes,'int16');
im_temp = zeros(N,M);
disp('Loading file....')
for loadframe = 1:numframes
    im(:,:,loadframe) = imread([tiff_path tiff_filename],'tiff',loadframe);
    disp(['Frame ' num2str(loadframe) '/' num2str(numframes)]);
end
disp('Done')
handles.im = im;
M = size(handles.im,1);
N = size(handles.im,2);
P = size(handles.im,3);
handles.im = reshape(handles.im,M*N,P);
disp('CURRENTLY LOADING IN AS DOUBLE')
handles.im = double(handles.im);
%%%%

axes(handles.guiAxes);
cla(handles.guiAxes);
imagesc(reshape(handles.im(:,1),M,N),'CDataMapping','scaled');
xlim([0 N]);
ylim([0 M]);
colormap(gray)
[imgMin imgMax] = caxis;
set(handles.imgMin,'string',imgMin);
set(handles.imgMax,'string',imgMax);
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
AP_cellROI_getTrace
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
    %set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on slider movement.
function imgSlider_Callback(hObject, eventdata, handles)
% hObject    handle to imgSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get frame number from slider, round appropriately
tiff_fullfilename = get(handles.tiffFile,'string');
frame = get(handles.imgSlider,'Value');
frame = round(frame);
set(handles.imgSlider,'Value',frame);

[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

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

% store ROIs in roilist edit box
h_autosort = figure;
AP_cellROIAutosort;
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

% make sure positions are updated % this should already be true
polywait = waitbar(0,'Saving traces');
for num_polygon = 1:length(polygon.ROI)
    %polygon.ROI{num_polygon} = getPosition(polygon.handles{num_polygon});
    roi_mask{num_polygon} = poly2mask(polygon.ROI{num_polygon}(:,1),polygon.ROI{num_polygon}(:,2),N,M);
    listselect = num_polygon;
    % get traces for all rois
    waitbar(num_polygon/length(polygon.ROI),polywait,'Saving ROI regions');
end

close(polywait);
disp('Not saving trace for now, uncomment line 766 if you want this')
%AP_cellROI_saveTrace
save([roi_savepath roi_savefile], 'roi_mask','polygon')%'roi_trace',


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
    polygon.handles{n_polygon} = line(polygon.ROI{n_polygon}(:,1),polygon.ROI{n_polygon}(:,2));
    set(polygon.handles{n_polygon}, 'ButtonDownFcn', {@selectROI, handles});
end

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
    delete(polygon.handles{oldROI})
    polygon.handles{oldROI} = line(polygon.ROI{oldROI}(:,1), polygon.ROI{oldROI}(:,2));
    set(polygon.handles{oldROI}, 'ButtonDownFcn', {@selectROI, handles});
end


% redraw all polygons in blue
currpoly = findobj(handles.guiAxes,'type','line');

for n_polygon = 1:length(polygon.ROI)
    set(currpoly(n_polygon), 'color', 'b');
end

% draw selected polygon in blue
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

% get current file and file to be compared
tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_orig_file tiff_orig_path] = uigetfile('*.tif');

im_curr(:,:)=double(imread(tiff_fullfilename,'tiff',1));
im_orig(:,:)=double(imread([tiff_orig_path tiff_orig_file],'tiff',1));

% DON'T HARDCODE THIS IN THE FUTURE
disp 'xy_ratio hardcoded as 4'
xy_ratio = 4;

% get estimate for dx dy
cc = normxcorr2(im_orig,im_curr);
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));

% maximize correlation through transformation
starting = [xpeak ypeak 0]; % dx dy theta
options = optimset('Algorithm','interior-point');
lb = [-100 -100/xy_ratio -pi/2];
ub = [100 100/xy_ratio pi/2];
% estimate transform of orig -> curr
display 'Finding best fit...'
estimates = fmincon(@(x) AP_affineAlign(x,im_curr,im_orig,xy_ratio),starting,[],[],[],[],lb,ub,[],options);
display 'Done'
% get estimates from best guess
dx = estimates(1);
dy = estimates(2);
theta = estimates(3);
[N M] = size(im_orig);

axes(handles.guiAxes);

% center image and constrain output to image size
[N M] = size(im_orig);

udata = [1 N] - (N+1)/2;% input image x
vdata = [1 M] - (M+1)/2;% input image y

xdata = udata;% output image x
ydata = vdata;% output image y

% define affine transform matrix
A = ...
    [cos(theta) sin(theta) 0;
    -sin(theta) cos(theta) 0;
    dx          dy          1];

tform_roi = maketform('affine',A);

% transform polygons according to image match
for i = 1:length(polygon.handles)
    delete(polygon.handles{i})
    % transform polygon, fix for xy_ratio
    polygon.ROI{i} = tformfwd(tform_roi,[polygon.ROI{i}(:,1) polygon.ROI{i}(:,2).*xy_ratio]);%, ...
    % 'UData',udata,'VData',vdata,'XData',xdata,'YData',ydata);
    % correct final ROIs for xy_ratio
    polygon.ROI{i}(:,2) = polygon.ROI{i}(:,2)./xy_ratio;

    % update polygon lines
    polygon.handles{i} = line(polygon.ROI{i}(:,1), polygon.ROI{i}(:,2));
    set(polygon.handles{i}, 'ButtonDownFcn', {@selectROI, handles});
end

% update roi list
set(handles.roiList,'UserData',polygon);


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
            set(handles.roi_numbers_h(currpoly),'Visible','on');
        end
    end
else
    for currpoly = 1:length(polygon.handles);
        set(polygon.handles{currpoly},'Visible','off');
        if isfield(handles,'roi_numbers_h')
            set(handles.roi_numbers_h(currpoly),'Visible','off');
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
im_orig = double(imread([orig_path orig_file],'tiff'));

% scale images to be 0-255 uint 8
im_orig = uint8((im_orig./(max(im_orig(:))))*255);
im_align = reshape(mean(handles.im,2),N,M);
im_align = uint8((im_align./(max(im_align(:))))*255);

[input_points base_points] = cpselect(im_orig,im_align,'Wait',true);
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
    [imgMin imgMax] = caxis;
    set(handles.imgMin,'string',imgMin);
    set(handles.imgMax,'string',imgMax);
    set(gca,'Children',flipud(get(gca,'Children')));

elseif stdImg == 1;
    if strcmp(class(handles.im),'int16')
        disp('Cannot do std at the moment - requires double');
    elseif strcmp(class(handles.im),'double')
        set(handles.imgSlider,'Enable','Off');
        axes(handles.guiAxes);
        im_h = findobj('type','image');
        delete(im_h(1));
        hold on;
        imagesc(reshape(nanstd(handles.im,0,2),N,M),'CDataMapping','scaled');
        colormap(gray)
        [imgMin imgMax] = caxis;
        set(handles.imgMin,'string',imgMin);
        set(handles.imgMax,'string',imgMax);
        %set(gca,'Children',flipud(get(gca,'Children')));
    end
end


% --- Executes on button press in helpROI.
function helpROI_Callback(hObject, eventdata, handles)
% hObject    handle to helpROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if ~isfield(handles,'im_smooth')
%     disp 'NOT SMOOTHING AT THE MOMENT'
%     M = size(handles.im,1);
%     N = size(handles.im,2);
%     P = size(handles.im,3);
%     im_long = reshape(handles.im,M*N,P);
%     im_smooth = im_long;
% %     im_smooth = zeros(size(handles.im));
% %     M = size(handles.im,1);
% %     N = size(handles.im,2);
% %     P = size(handles.im,3);
% %     im_long = reshape(handles.im,M*N,P);
% %     im_smooth = zeros(M*N,P);
% %     disp('Smoothing movie...')
% %     parfor i = 1:M*N;
% %         im_smooth(i,:) = smooth(im_long(i,:),10);
% %     end
% %     disp('done')
%     handles.im_smooth = im_smooth;
%     % update handles
%     guidata(hObject, handles);
% end
std_thresh = str2num(get(handles.std_thresh,'string'));
tiff_fullfilename = get(handles.tiffFile,'string');
[tiff_path tiff_filename] = fileparts(tiff_fullfilename);
imageinfo=imfinfo([tiff_path filesep tiff_filename],'tiff');
M=imageinfo(1).Width;
N=imageinfo(1).Height;

axes(handles.guiAxes);
[all_x all_y] = ginput(1);
for cell_num = 1:length(all_x);
    x = round(all_x(cell_num));
    y = round(all_y(cell_num));
    ind = sub2ind([N M],y,x);
    corrcoef_grid = [];
    
    if strcmp(class(handles.im),'double')
        corrcoef_grid = [handles.im(ind,:)*handles.im']' ./ sum(handles.im,2);
    elseif strcmp(class(handles.im),'int16')
        % do matrix math in a for loop :(
        tic
        for num_pixel = 1:size(handles.im,1)
            corrcoef_grid(num_pixel) = int16(sum(double(handles.im(ind,:)).*double(handles.im(num_pixel,:))) / double(sum(handles.im(num_pixel,:))));
            disp(['Finding dot product ' num2str(num_pixel) '/' num2str(size(handles.im,1))]);
        end
        toc
    end
    %

    corrcoef_grid_reshape = reshape(corrcoef_grid,N,M);
    % threshold the corrcoef_grid_reshape
    corrcoef_threshold = int16(nanmean(corrcoef_grid) + std_thresh*nanstd(corrcoef_grid));
    binaryImg = [];
    binaryImg = corrcoef_grid_reshape >= corrcoef_threshold; % threshold image
    if any(any(binaryImg))
        binaryCell = [];
        binaryCell = bwselect(binaryImg,x,y); % select cell that was clicked
        if any(any(binaryCell))
            first_nonzero = find(binaryCell > 0,1);
            [y_nonzero x_nonzero] = ind2sub([N M],first_nonzero);
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
end
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
    %set(gca,'Children',flipud(get(gca,'Children')));
end


% --- Executes on button press in showBehavior.
function showBehavior_Callback(hObject, eventdata, handles)
% hObject    handle to showBehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showBehavior

n_channels = 1; % true at the moment for all motion corrected files
n_frames = 1000; %get this later

showBehavior = get(handles.showBehavior,'value');

if showBehavior == 1
    [bhv_filename bhv_path] = uigetfile('*','Behavior file');
    bhv = load([bhv_path bhv_filename],'-MAT');

    saved_history = [];
    saved_history = bhv.saved_history; % this in for now, somewhere below fix

    [acq_filename_all acq_path] = uigetfile('*.xsg','Choose XSG files','Multiselect','on');

    disp 'Finding events...'
    
    left_down_all = [];
    left_up_all = [];
    reward_all = [];
    unrewarded_all = [];
    cue_onset_all = [];
    frames_left_down_all = [];
    frames_left_up_all = [];
    frames_reward_all = [];
    frames_unrewarded_all = [];
    frames_cue_onset_all = [];
    timeDiff = [];
    reward_trial = [];
    unrewarded_trial = [];

    for i = 1:length(acq_filename_all);

        % get current acq filename
        acq_filename = [acq_path acq_filename_all{i}];

        % this is hard coded until the header problem is fixed
        % % %     % pull out time for each for each frame from header
        % % %     img_info = imfinfo([img_pathname img_filename]);
        % % %     msPerLine_start = strfind(img_info(1).ImageDescription, 'state.acq.msPerLine')+20;
        % % %     msPerLine_end = strfind(img_info(1).ImageDescription, 'state.acq.fillFraction')-2;
        % % %     msPerLine = str2num(img_info(1).ImageDescription(msPerLine_start:msPerLine_end));
        % % %     lines = img_info(1).Height;
        % % %     fpms = msPerLine*lines;

        msPerLine = 1.24;
        lines = 128;
        fpms = msPerLine*lines;

        % Get trials in seconds since started
        TrialNumList = TK_ReadBitCode([acq_filename]);
        close(gcf)
        % put trial time in ms
        TrialNumList(2,:) = TrialNumList(2,:)*1000;

        left_down = [];
        left_up = [];
        reward = [];
        unrewarded = [];
        cue_onset = [];

        for trial = TrialNumList(1,:)
            if trial > length(bhv.saved_history.RewardsSection_LastTrialEvents);
                continue
            end
            bhv_current = bhv.saved_history.RewardsSection_LastTrialEvents{trial,1};
            % fixes possible bug of putting in old timepoints at the end
            LastGoodTrial = find(diff(bhv_current(:,3)) < 0,1);
            if ~isempty(LastGoodTrial)
                bhv_current = bhv_current(1:LastGoodTrial,:);
            end
            % get time lag of behavior to acq (t_bhv + timeDiff = t_acq) in ms
            trialStart_bhv = bhv_current(find(bhv_current(:,1) == 101 & bhv_current(:,2) == 0,1,'last'),3)*1000;
            lag_time = TrialNumList(2,find(TrialNumList(1,:) == trial)) - trialStart_bhv;
            timeDiff(trial,1) = lag_time;
            % make a second column for which acquisition you're in to frame shift
            timeDiff(trial,2) = i*ones(length(lag_time),1);
            % time for all lever up and down (left down = 5, left up = 6) in ms
            left_down = [left_down; bhv_current(find(bhv_current(:,2) == 5),3)*1000 + timeDiff(trial)];
            left_up = [left_up; bhv_current(find(bhv_current(:,2) == 6),3)*1000 + timeDiff(trial)];
            % time for all rewarded press (state 44) and unrewarded press (state 45)
            reward = [reward; bhv_current(find(bhv_current(:,1) == 44,1),3)*1000 + timeDiff(trial,1)];
            unrewarded = [unrewarded; bhv_current(find(bhv_current(:,1) == 45,1),3)*1000 + timeDiff(trial,1)];
            % time for cue onset (state 41)
            cue_onset = [cue_onset; bhv_current(find(bhv_current(:,1) == 41,1),3)*1000 + timeDiff(trial,1)];
        end

        % concatenate all events, multiply by img number to get total time

        left_down_all = [left_down_all; left_down];
        left_up_all = [left_up_all; left_up];
        unrewarded_all = [unrewarded_all; unrewarded];
        reward_all = [reward_all; reward];
        cue_onset_all = [cue_onset_all; cue_onset];

        frames_left_down = ceil(left_down(~isnan(left_down))/fpms)';
        frames_left_up = ceil(left_up(~isnan(left_up))/fpms)';
        frames_reward = ceil(reward(~isnan(reward))/fpms)';
        frames_unrewarded = ceil(unrewarded(~isnan(unrewarded))/fpms)';
        frames_cue_onset = ceil(cue_onset(~isnan(cue_onset))/fpms)';

        % round up for number of channels
        frames_left_down = (1-mod(frames_left_down,n_channels))+frames_left_down;
        frames_left_up = (1-mod(frames_left_up,n_channels))+frames_left_up;
        frames_cue_onset = (1-mod(frames_cue_onset,n_channels))+frames_cue_onset;
        frames_reward = (1-mod(frames_reward,n_channels))+frames_reward;
        frame_unrewarded = (1-mod(frames_unrewarded,n_channels))+frames_unrewarded;

        % Get frame for each event
        frames_left_down_all = [frames_left_down_all frames_left_down+n_frames*(i-1)];
        frames_left_up_all = [frames_left_up_all frames_left_up+n_frames*(i-1)];
        frames_reward_all = [frames_reward_all frames_reward+n_frames*(i-1)];
        frames_unrewarded_all = [frames_unrewarded_all frames_unrewarded+n_frames*(i-1)];
        frames_cue_onset_all = [frames_cue_onset_all frames_cue_onset+n_frames*(i-1)];

    end

    cd ..
    
    axes(handles.behaviorPlot)
    hold on
    plot(frames_left_down_all,1,'.b')
    plot(frames_left_up_all,2,'.k')
    plot(frames_reward_all,3,'.g')
    plot(frames_unrewarded_all,3,'.r')
    plot(frames_cue_onset_all,4,'.c')
    ylim([1 4])
    set(handles.behaviorPlot,'YTickLabel',{'Down','Up','Reward','Cue'})
    disp 'Done'
elseif showBehavior == 0;
    axes(handles.behaviorPlot)
    clear axes
end



function std_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to std_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std_thresh as text
%        str2double(get(hObject,'String')) returns contents of std_thresh as a double


% --- Executes during object creation, after setting all properties.
function std_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
