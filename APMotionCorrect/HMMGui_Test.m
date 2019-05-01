function varargout = HMMGui_Test(varargin)
% HMMGUI_TEST MATLAB code for HMMGui_Test.fig
%      HMMGUI_TEST, by itself, creates a new HMMGUI_TEST or raises the existing
%      singleton*.
%
%      H = HMMGUI_TEST returns the handle to a new HMMGUI_TEST or the handle to
%      the existing singleton*.
%
%      HMMGUI_TEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HMMGUI_TEST.M with the given input arguments.
%
%      HMMGUI_TEST('Property','Value',...) creates a new HMMGUI_TEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HMMGui_Test_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HMMGui_Test_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HMMGui_Test

% Last Modified by GUIDE v2.5 28-Jan-2011 17:25:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HMMGui_Test_OpeningFcn, ...
    'gui_OutputFcn',  @HMMGui_Test_OutputFcn, ...
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


% --- Executes just before HMMGui_Test is made visible.
function HMMGui_Test_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HMMGui_Test (see VARARGIN)

% Choose default command line output for HMMGui_Test
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HMMGui_Test wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HMMGui_Test_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function MaxDx_Callback(hObject, eventdata, handles)
% hObject    handle to MaxDx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxDx as text
%        str2double(get(hObject,'String')) returns contents of MaxDx as a double


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

% Hints: get(hObject,'String') returns contents of MaxDy as text
%        str2double(get(hObject,'String')) returns contents of MaxDy as a double


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

% Hints: get(hObject,'String') returns contents of Lambda as text
%        str2double(get(hObject,'String')) returns contents of Lambda as a double


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




function savedir_Callback(hObject, eventdata, handles)
% hObject    handle to savedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savedir as text
%        str2double(get(hObject,'String')) returns contents of savedir as a double


% --- Executes during object creation, after setting all properties.
function savedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savedir (see GCBO)
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
[path]=uigetdir;
set(handles.savedir,'string',path);


function TIFFFile_Callback(hObject, eventdata, handles)
% hObject    handle to TIFFFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TIFFFile as text
%        str2double(get(hObject,'String')) returns contents of TIFFFile as a double


% --- Executes during object creation, after setting all properties.
function TIFFFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TIFFFile (see GCBO)
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
[filename,path]=uigetfile('*.*','Pick file to be corrected','Multiselect','on');
set(handles.TIFFFile,'string',[path filename]);


function TIFFRef_Callback(hObject, eventdata, handles)
% hObject    handle to TIFFRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TIFFRef as text
%        str2double(get(hObject,'String')) returns contents of TIFFRef as a double


% --- Executes during object creation, after setting all properties.
function TIFFRef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TIFFRef (see GCBO)
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
[filename,path]=uigetfile('*.*','Pick reference file');
set(handles.TIFFRef,'string',[path filename]);

% --- Executes on button press in HMM.
function HMM_Callback(hObject, eventdata, handles)
% hObject    handle to HMM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
filelist = get(handles.filelist,'UserData');
if exist('filelist')
    
    % loop through all files in list
    for i = 1:size(filelist,1)
        
        % Get paths/filenames from GUI
        savedir = filelist{i}{1,1};
        filename = filelist{i}{2,1};
        reffilename = filelist{i}{3,1};
        [filename_path,filename_file,filename_ext] = fileparts(filename);
        
        %filename = get(handles.TIFFFile,'string');
        %savedir = get(handles.savedir,'string');
        %reffilename = get(handles.TIFFRef,'string');
        
        % Get settings from GUI
        autoRef = get(handles.AutomaticRef,'value');
        autoLambda = get(handles.AutomaticLambda,'value');
        playback =  get(handles.Playback,'value');
        maxdx = str2num(get(handles.MaxDx,'string'));
        maxdy = str2num(get(handles.MaxDy,'string'));
        channels = get(handles.channels,'value');
        turboreg = get(handles.turboreg,'value');
        lambda = str2num(get(handles.Lambda,'string'));

        
        % make filename for future saved HMM files
        savecorrected = [savedir filesep filename_file '_HMMAPTK_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
        
        % Process with chosen reference picture, or autopick one every 30s
        if autoRef == 1;
            split_createPI(filename,savedir,maxdx,maxdy,channels);
            load([savedir filesep filename_file '_PI.mat']);
        elseif autoRef == 0;
            split_createPI_ManRef(filename,savedir,reffilename,0,maxdx,maxdy,channels,lambda);
            load([savedir filesep filename_file '_Lambda_' num2str(lambda) '_PI.mat']);
        end
        
        % Calculate lambda or use given lambda
        if autoLambda == 0;
            lambda = str2num(get(handles.Lambda,'string'));
            markov_on_PIsave
        elseif autoLambda == 1;
            maxlambda;
        end
        
        % Play movie back after processing
        %if playback == 1;
            [fixeddata,countdata]=tk_playback_markov2(filename,offsets,edgebuffer,playback,lambda);
        %end
        %if playback ~= 1;
        %    [fixeddata,countdata]=playback_markov(chone,offsets,edgebuffer,0);
        %end    
        timeHMM = toc;
        im_c = permute(fixeddata,[2 3 1]);
        save(savecorrected, 'im_c','timeHMM')
    end
end

% --- Executes on button press in AutomaticRef.
function AutomaticRef_Callback(hObject, eventdata, handles)
% hObject    handle to AutomaticRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutomaticRef
valueAutomaticRef = get(handles.AutomaticRef,'value');
if valueAutomaticRef == 0
    set(handles.TIFFRef,'Enable','on')
    set(handles.RefBrowse,'Enable','on')
else
    set(handles.TIFFRef,'Enable','off')
    set(handles.RefBrowse,'Enable','off')
end

% --- Executes on button press in AutomaticLambda.
function AutomaticLambda_Callback(hObject, eventdata, handles)
% hObject    handle to AutomaticLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AutomaticLambda
valueAutomaticLambda = get(handles.AutomaticLambda,'value');
if valueAutomaticLambda == 0
    set(handles.Lambda,'Enable','on')
else
    set(handles.Lambda,'Enable','off')
end


% --- Executes on button press in Playback.
function Playback_Callback(hObject, eventdata, handles)
% hObject    handle to Playback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Playback


% --- Executes on selection change in channels.
function channels_Callback(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns channels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from channels


% --- Executes during object creation, after setting all properties.
function channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in turboreg.
function turboreg_Callback(hObject, eventdata, handles)
% hObject    handle to turboreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of turboreg


% --- Executes on selection change in filelist.
function filelist_Callback(hObject, eventdata, handles)
% hObject    handle to filelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filelist


% --- Executes during object creation, after setting all properties.
function filelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addfile.
function addfile_Callback(hObject, eventdata, handles)
% hObject    handle to addfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = get(handles.TIFFFile,'string');
savedir = get(handles.savedir,'string');
reffilename = get(handles.TIFFRef,'string');
%WARNING
disp 'this way of calling multiple files probz doesnt work - check it out'
if iscell(filename);
    for filenum = 2:length(filename);
        [filepath filebase ext] = fileparts(filename{filenum});
        filepath = filename{1};
        % update listbox
        currentlist = get(handles.filelist,'string');
        currentud = get(handles.filelist,'UserData');
        set(handles.filelist,'string',strvcat(currentlist,filebase));
        % update UserData
        newRow = size(currentud,1)+1;
        currentud{newRow,1} = {savedir;[filepath filebase];reffilename};
        set(handles.filelist,'UserData',currentud)
    end
else
    [filepath filebase ext] = fileparts(filename);
    % update listbox
    currentlist = get(handles.filelist,'string');
    currentud = get(handles.filelist,'UserData');
    set(handles.filelist,'string',strvcat(currentlist,filebase));
    % update UserData
    newRow = size(currentud,1)+1;
    currentud{newRow,1} = {savedir;filename;reffilename};
    set(handles.filelist,'UserData',currentud)
end


% --- Executes on button press in deletefile.
function deletefile_Callback(hObject, eventdata, handles)
% hObject    handle to deletefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listselect = get(handles.filelist,'value');
% delete from listbox
currentlist = get(handles.filelist,'string');
currentlist(listselect,:) = [];
% delete from usrdata
currentud = get(handles.filelist,'UserData');
currentud(listselect,:) = [];

set(handles.filelist,'value',1);
set(handles.filelist,'string',currentlist);
set(handles.filelist,'UserData',currentud);


% --- Executes on button press in jk_mc.
function jk_mc_Callback(hObject, eventdata, handles)
% hObject    handle to jk_mc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load in data

tic

filelist = get(handles.filelist,'UserData');
if exist('filelist')
    
    % loop through all files in list
    for i = 1:size(filelist,1)
        
        % Get paths/filenames from GUI
        savedir = filelist{i}{1,1};
        filename = filelist{i}{2,1};
        reffilename = filelist{i}{3,1};
        [filename_path, filename_file, filename_ext] = fileparts(filename);
        
        % create filename for future saved Kerr files
        savename = [savedir filesep filename_file '_KerrAP'];
      
        % Get settings from GUI
        autoRef = get(handles.AutomaticRef,'value');
        autoLambda = get(handles.AutomaticLambda,'value');
        playback =  get(handles.Playback,'value');
        maxdx = str2num(get(handles.MaxDx,'string'));
        maxdy = str2num(get(handles.MaxDy,'string'));
        channels = get(handles.channels,'value');
        turboreg = get(handles.turboreg,'value');
        
        % This loads files assuming that they're 1 channel
        disp 'This is assuming 1 channel and reference = ref at the moment as .mat'
        
        imageinfo=imfinfo(filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;    
        im_s=zeros(N,M,numframes);
        
        % Read in tiff file to be corrected
        for u = 1:numframes
        im_s(:,:,u)=imread(filename,'tiff',u);
        end
        
        % If you've got NaN, make 0 (keep in mind for correlation)
        im_s(im_s<1)=0;
        load(reffilename);
        im_t = ref;
        %im_t = mean(im_s,3);
        
        % Do Kerr correction ONLY FOR 'framelimit' FRAMES AT A TIME
        framelimit = 600; % frames allowed in each processing chunk
        disp 'doing framelimit of 600'
      
        numframes > framelimit;
        snips = ceil(numframes/framelimit);

        for g = 1:snips;
                
            if g ~= snips
                im_s_snips = im_s(:,:,(g-1)*framelimit+1:(g-1)*framelimit+framelimit);
            else
                im_s_snips = im_s(:,:,(g-1)*framelimit+1:end);
            end
            
            [im_c dx_r dy_r E] = imreg_greenberg(im_s_snips, im_t, 1); 
                        
            timeK = toc;   
            savename_current = [savename '_' num2str(g) '-' num2str(snips) '.mat'];
            save(savename_current,'im_c','timeK');
            clear im_s_snips im_c         
        end
    end
end


% --- Executes on button press in turbreg.
function turbreg_Callback(hObject, eventdata, handles)
% hObject    handle to turbreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filelist = get(handles.filelist,'UserData');
if exist('filelist')
    
    % loop through all files in list
    for i = 1:size(filelist,1)
        
        % Get paths/filenames from GUI
        savedir = filelist{i}{1,1};
        filename = filelist{i}{2,1};
        reffilename = filelist{i}{3,1};
        
        [sourcepath, sourcename, sourceext] = fileparts(filename);
        
        % create filename for future saved Turboreg files
        savename = [savedir filesep sourcename '_TurboRegAP'];
        
        % Get settings from GUI
        autoRef = get(handles.AutomaticRef,'value');
        autoLambda = get(handles.AutomaticLambda,'value');
        playback =  get(handles.Playback,'value');
        maxdx = str2num(get(handles.MaxDx,'string'));
        maxdy = str2num(get(handles.MaxDy,'string'));
        channels = get(handles.channels,'value');
        turboreg = get(handles.turboreg,'value');
        
        % This loads files assuming that they're 1 channel
        disp 'This is assuming 1 channel and reference = ref at the moment as .mat'
        
        imageinfo=imfinfo(filename,'tiff');
        numframes=length(imageinfo);
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        im_s=zeros(N,M,numframes);
        
        % Read in tiff file to be corrected
        for u = 1:numframes
            im_s(:,:,u)=imread(filename,'tiff',u);
        end
        
        
        % If you've got NaN, make 0 (keep in mind for correlation)
        im_s(im_s<1)=0;
        load(reffilename);
        im_t = ref;
        %im_t = mean(im_s,3);
        
        
        javaaddpath 'C:\Matlab2010B\MATHWORKS_R2010B\java\jar\ij.jar'
        javaaddpath 'C:\Matlab2010B\MATHWORKS_R2010B\java\jar\mij.jar'
        javaaddpath 'C:\Matlab2010B\MATHWORKS_R2010B\java'
        
        % overwrite old file with new correction
        
        % Keep in mind that you need a averaged image called ref
        disp 'Modified this to use ref avgd before hand not use autoref VARIABLE MUST BE CALLED REF'
        if autoRef == 0;
            
            disp 'if not tiff, assumed variable called im_c'
            Turboreg_AP(reffilename, filename, savename);%(targetfilename, sourcefilename, savefilename)
        end
    end
end

