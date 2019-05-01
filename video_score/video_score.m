function varargout = video_score(varargin)
% VIDEO_SCORE MATLAB code for video_score.fig
%      VIDEO_SCORE, by itself, creates a new VIDEO_SCORE or raises the existing
%      singleton*.
%
%      H = VIDEO_SCORE returns the handle to a new VIDEO_SCORE or the handle to
%      the existing singleton*.
%
%      VIDEO_SCORE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEO_SCORE.M with the given input arguments.
%
%      VIDEO_SCORE('Property','Value',...) creates a new VIDEO_SCORE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before video_score_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to video_score_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help video_score

% Last Modified by GUIDE v2.5 23-Nov-2013 01:23:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @video_score_OpeningFcn, ...
                   'gui_OutputFcn',  @video_score_OutputFcn, ...
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


% --- Executes just before video_score is made visible.
function video_score_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to video_score (see VARARGIN)

% Choose default command line output for video_score
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes video_score wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = video_score_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function keyPress(currentObject, eventdata, hObject)
% use shortcuts when image is selected
handles = guidata(hObject);
switch eventdata.Key
    case 'a'
        if length(handles.light_on) < length(handles.light_off)+1
            handles.light_on = [handles.light_on toc(handles.master_clock)];
            disp('Light on')
        end
    case 'l'
        if length(handles.movement_on) < length(handles.movement_off)+1
            handles.movement_on = [handles.movement_on toc(handles.master_clock)];
            disp('Movement on')
        end
    case 'space'
        if ~handles.clock_running
            handles.actx.controls.play;
            handles.clock_running = true;
            handles.master_clock = tic;
            disp('Clock started');
        else
            handles.actx.controls.stop;
            handles.clock_running = false;
            handles.total_time = toc(handles.master_clock);
            disp('Clock stopped');
            close(handles.movie_fig)
        end
end
guidata(hObject,handles);

function keyRelease(currentObject, eventdata, hObject)
% use shortcuts when image is selected
handles = guidata(hObject);
switch eventdata.Key
    case 'a'
        handles.light_off = [handles.light_off toc(handles.master_clock)];
        disp('Light off')
    case 'l'
        handles.movement_off = [handles.movement_off toc(handles.master_clock)];
        disp('Movement off')
end
guidata(hObject,handles);


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
times.master_clock = handles.master_clock;
times.total_time = handles.total_time;

times.light_on = handles.light_on;
times.light_off = handles.light_off;

times.movement_on = handles.movement_on;
times.movement_off = handles.movement_off;

user_initials = input('Initials ? ','s');
savename = [handles.movie_filename(1:end-4) ...
    '_score_' ...
    user_initials];

[savefile savepath] = ...
    uiputfile(savename);
save([savepath savefile],'times');


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile('*.avi','Select movie');
handles.movie_filename = filename;
handles.movie_fig = figure;

% set button presses for figure
set(gcf,'KeyPressFcn', {@keyPress, hObject}, ...
    'KeyReleaseFcn', {@keyRelease, hObject}); 

% get and load movie
handles.actx = actxcontrol('WMPlayer.OCX.7', [0 0 500 500],handles.movie_fig);
set(handles.actx.settings,'autoStart',0)
handles.actx.URL=[pathname filename];

handles.light_on = [];
handles.light_off = [];

handles.movement_on = [];
handles.movement_off = [];

handles.clock_running = false;

guidata(hObject, handles);

