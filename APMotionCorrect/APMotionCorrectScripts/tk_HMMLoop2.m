function varargout = tk_HMMLoop2(varargin)
% TK_HMMLOOP2 M-file for tk_HMMLoop2.fig
%      TK_HMMLOOP2, by itself, creates a new TK_HMMLOOP2 or raises the existing
%      singleton*.
%
%      H = TK_HMMLOOP2 returns the handle to a new TK_HMMLOOP2 or the handle to
%      the existing singleton*.
%
%      TK_HMMLOOP2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TK_HMMLOOP2.M with the given input arguments.
%
%      TK_HMMLOOP2('Property','Value',...) creates a new TK_HMMLOOP2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tk_HMMLoop2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tk_HMMLoop2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tk_HMMLoop2

% Last Modified by GUIDE v2.5 31-Dec-2010 03:23:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tk_HMMLoop2_OpeningFcn, ...
                   'gui_OutputFcn',  @tk_HMMLoop2_OutputFcn, ...
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


% --- Executes just before tk_HMMLoop2 is made visible.
function tk_HMMLoop2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tk_HMMLoop2 (see VARARGIN)
global setting
% Choose default command line output for tk_HMMLoop2
handles.output = hObject;

setting.maxdx = 20;
set(handles.MaxDx,'String',setting.maxdx)
setting.maxdy = 5;
set(handles.MaxDy,'String',setting.maxdy)
setting.lambda = 0.05;
set(handles.Lambda,'String',setting.lambda)
set(handles.LoadMode0,'Value',1)
setting.numchannels = 1;
set(handles.greenonly,'Value',1)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tk_HMMLoop2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tk_HMMLoop2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadMode0.
function LoadMode0_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode0
if get(hObject,'Value')==1
    setting.loadmode = 0;
    set(handles.LoadMode1, 'Value', 0)
    set(handles.LoadMode2, 'Value', 0)
    set(handles.LoadMode3, 'Value', 0)
    set(handles.LoadMode4, 'Value', 0)
    set(handles.LoadMode5, 'Value', 0)
end

% --- Executes on button press in LoadMode1.
function LoadMode1_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode1
if get(hObject,'Value')==1
    setting.loadmode = 1;
    set(handles.LoadMode0, 'Value', 0)
    set(handles.LoadMode2, 'Value', 0)
    set(handles.LoadMode3, 'Value', 0)
    set(handles.LoadMode4, 'Value', 0)
    set(handles.LoadMode5, 'Value', 0)
end

% --- Executes on button press in LoadMode2.
function LoadMode2_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode2
if get(hObject,'Value')==1
    setting.loadmode = 2;
    set(handles.LoadMode0, 'Value', 0)
    set(handles.LoadMode1, 'Value', 0)
    set(handles.LoadMode3, 'Value', 0)
    set(handles.LoadMode4, 'Value', 0)
    set(handles.LoadMode5, 'Value', 0)
end

% --- Executes on button press in LoadMode3.
function LoadMode3_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode3
if get(hObject,'Value')==1
    setting.loadmode = 3;
    set(handles.LoadMode0, 'Value', 0)
    set(handles.LoadMode1, 'Value', 0)
    set(handles.LoadMode2, 'Value', 0)
    set(handles.LoadMode4, 'Value', 0)
    set(handles.LoadMode5, 'Value', 0)
end

% --- Executes on button press in LoadMode4.
function LoadMode4_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode4
if get(hObject,'Value')==1
    setting.loadmode = 4;
    set(handles.LoadMode0, 'Value', 0)
    set(handles.LoadMode1, 'Value', 0)
    set(handles.LoadMode2, 'Value', 0)
    set(handles.LoadMode3, 'Value', 0)
    set(handles.LoadMode5, 'Value', 0)
end



% --- Executes on button press in LoadMode5.
function LoadMode5_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMode5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of LoadMode5
if get(hObject,'Value')==1
    setting.loadmode = 5;
    set(handles.LoadMode0, 'Value', 0)
    set(handles.LoadMode1, 'Value', 0)
    set(handles.LoadMode2, 'Value', 0)
    set(handles.LoadMode3, 'Value', 0)
    set(handles.LoadMode4, 'Value', 0)
end


function MaxDx_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hints: get(hObject,'String') returns contents of MaxDx as text
%        str2double(get(hObject,'String')) returns contents of MaxDx as a double
setting.maxdx = str2double(get(hObject,'String'));



% --- Executes during object creation, after setting all properties.
function MaxDx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxDy_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hints: get(hObject,'String') returns contents of MaxDy as text
%        str2double(get(hObject,'String')) returns contents of MaxDy as a double
setting.maxdy = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function MaxDy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxDy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Lambda_Callback(hObject, eventdata, handles)
% hObject    handle to Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hints: get(hObject,'String') returns contents of Lambda as text
%        str2double(get(hObject,'String')) returns contents of Lambda as a double
setting.lambda = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function Lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in greenonly.
function greenonly_Callback(hObject, eventdata, handles)
% hObject    handle to greenonly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of greenonly
if get(hObject,'Value')==1
    setting.numchannels = 1;
    set(handles.twochannel, 'Value', 0)
end

% --- Executes on button press in twochannel.
function twochannel_Callback(hObject, eventdata, handles)
% hObject    handle to twochannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% Hint: get(hObject,'Value') returns toggle state of twochannel
if get(hObject,'Value')==1
    setting.numchannels = 2;
    set(handles.greenonly, 'Value', 0)
end



function DirName_Callback(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirName as text
%        str2double(get(hObject,'String')) returns contents of DirName as a double


% --- Executes during object creation, after setting all properties.
function DirName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in DirBrowse.
function DirBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to DirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
% setting.filedir  = uigetdir('C:\Program Files\MATLAB\R2007b\work\HMM','Pick a source directory');
setting.filedir  = uigetdir('E:\ImagingData','Pick a source directory');
set(handles.DirName,'String',setting.filedir)



function FileName_Callback(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double


% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FileBrowse.
function FileBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to FileBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
%     [filename,pathname]=uigetfile('C:\Program Files\MATLAB\R2007b\work\HMM','pick your tif file','*.tif');
        [filename,pathname]=uigetfile('E:\ImagingData','pick your tif file','*.tif');
    setting.fullfilename=[pathname filesep filename];
set(handles.FileName,'String',setting.fullfilename)

function RefName_Callback(hObject, eventdata, handles)
% hObject    handle to RefName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RefName as text
%        str2double(get(hObject,'String')) returns contents of RefName as a double


% --- Executes during object creation, after setting all properties.
function RefName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RefName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RefBrowse.
function RefBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to RefBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
%     [filename,pathname]=uigetfile('C:\Program Files\MATLAB\R2007b\work\HMM','pick your tif file','*.tif');
    [filename,pathname]=uigetfile('E:\ImagingData','pick your tif file','*.tif');
    setting.fullreferencefilename=[pathname filesep filename];
set(handles.RefName,'String',setting.fullreferencefilename)


% --- Executes on button press in HMM.
function HMM_Callback(hObject, eventdata, handles)
% hObject    handle to HMM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
tk_HMM(setting.fullfilename,setting.fullreferencefilename,setting.loadmode,setting.maxdx,...
    setting.maxdy,setting.numchannels,setting.lambda)


% --- Executes on button press in HMMLoop.
function HMMLoop_Callback(hObject, eventdata, handles)
% hObject    handle to HMMLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global setting
tk_HMMLoop(setting.filedir,setting.fullreferencefilename,setting.loadmode,setting.maxdx,...
    setting.maxdy,setting.numchannels,setting.lambda)


% --- Executes on key press with focus on HMM and none of its controls.
function HMM_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to HMM (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
