function varargout = APMotionCorrect(varargin)
% APMOTIONCORRECT MATLAB code for APMotionCorrect.fig
%      APMOTIONCORRECT, by itself, creates a new APMOTIONCORRECT or raises the existing
%      singleton*.
%
%      H = APMOTIONCORRECT returns the handle to a new APMOTIONCORRECT or the handle to
%      the existing singleton*.
%
%      APMOTIONCORRECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APMOTIONCORRECT.M with the given input arguments.
%
%      APMOTIONCORRECT('Property','Value',...) creates a new APMOTIONCORRECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before APMotionCorrect_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to APMotionCorrect_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help APMotionCorrect

% Last Modified by GUIDE v2.5 04-May-2011 00:36:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @APMotionCorrect_OpeningFcn, ...
    'gui_OutputFcn',  @APMotionCorrect_OutputFcn, ...
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


% --- Executes just before APMotionCorrect is made visible.
function APMotionCorrect_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to APMotionCorrect (see VARARGIN)

% Choose default command line output for APMotionCorrect
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes APMotionCorrect wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = APMotionCorrect_OutputFcn(hObject, eventdata, handles)
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
[filename,path]=uigetfile('*.*','Pick reference file','Multiselect','on');
set(handles.TIFFRef,'string',[path filename]);


% --- Executes on button press in HMM.
function HMM_Callback(hObject, eventdata, handles)
% hObject    handle to HMM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get initial settings
master_filelist = get(handles.filelist,'UserData');
autoRef = get(handles.AutomaticRef,'value');
parsedfile = get(handles.parsedfile,'value');
channels = get(handles.channels,'value');
multisession = get(handles.multisession,'value');
autoLambda = get(handles.AutomaticLambda,'value');
playback =  get(handles.Playback,'value');
maxdx = str2num(get(handles.MaxDx,'string'));
maxdy = str2num(get(handles.MaxDy,'string'));
channels = get(handles.channels,'value');
turboreg = get(handles.turboreg,'value');
lambda = str2num(get(handles.Lambda,'string'));
lambda_string = get(handles.Lambda,'string');
% record only values after decimal to avoid filename errors
lambda_string = lambda_string(lambda_string ~= '.');
motionBias = get(handles.motionBias,'value');

if multisession == 1;
    num_session = max(master_filelist{end}{4});
elseif multisession == 0;
    filelist = master_filelist;
    num_session = 1;
end

for session = 1:num_session
    
    if multisession == 1;
        session_index = cellfun(@(x) x{4}==session,master_filelist);
        filelist = master_filelist(session_index);
    end
    
    % create parsed ref if selected
    if parsedfile == 1 && autoRef == 1;
        APMotionCorrect_parsedRef
    end
    
    if exist('filelist')
        
        % loop through files with turboreg if selected
        if turboreg == 1 ;
            
            APMotionCorrect_turboreg
            
            % create new parsed file avg if turboreg
            if parsedfile == 1 && autoRef == 1;
                APMotionCorrect_parsedTurboregRef
            end
            
        end
                
        % loop through files with HMM
        APMotionCorrect_HMM
        
    end
    
    % get rid of temporary ref files
    if exist([savedir filesep 'parsedAverage.tif'],'file')
        delete([savedir filesep 'parsedAverage.tif']);
    end
    if exist([savedir filesep 'parsedTurboregAverage.tif'],'file')
        delete([savedir filesep 'parsedTurboregAverage.tif']);
    end
    if exist([savedir filesep 'currfileAverage.tif'],'file')
        delete([savedir filesep 'currfileAverage.tif']);
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
multisession = get(handles.multisession,'value');

if ~iscell(filename);
    filename = mat2cell(filename);
end

for filenum = 1:length(filename);
    if length(filename) == 1;
        [filepath filebase ext] = fileparts(filename{filenum});
        filepath = [filepath filesep];
    else
        if filenum == 1;
            continue
        else
            [filepath filebase ext] = fileparts(filename{filenum});
            filepath = filename{1};
        end
    end
    % update listbox
    
    % do session numbers if multisession
    if multisession == 1
        currentlist = get(handles.filelist,'string');
        currentud = get(handles.filelist,'UserData');
        if length(filename) > 1
            if isempty(currentlist);
                n_group = 1;
            elseif ~isempty(currentlist) && filenum == 2;
                n_group = currentud{end}{end} + 1;
            elseif ~isempty(currentlist) && filenum > 2;
                n_group = currentud{end}{end};
            end
        end
        listboxContent = strvcat(currentlist,['SESSION ' num2str(n_group) ': ' filebase]);
        set(handles.filelist,'string',listboxContent);
        % update UserData
        newRow = size(currentud,1)+1;
        currentud{newRow,1} = {savedir;[filepath filebase];reffilename;n_group};
        set(handles.filelist,'UserData',currentud)
        
        % don't add session numbers if not multisession
    elseif multisession ==  0;
        currentlist = get(handles.filelist,'string');
        currentud = get(handles.filelist,'UserData');
        listboxContent = strvcat(currentlist,filebase);
        set(handles.filelist,'string',listboxContent);
        % update UserData
        newRow = size(currentud,1)+1;
        currentud{newRow,1} = {savedir;[filepath filebase];reffilename};
        set(handles.filelist,'UserData',currentud)
    end
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
        
        if autoRef == 0;
            im_t = imread(reffilename,'tiff');
        elseif autoRef == 1;
            im_t = mean(im_s,3);
        end
        
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
            h_waitbar = waitbar(0,'Saving corrected file');
            for frame = 1:length(im_c,3);
                a=uint16(round(im_c(:,:,frame)));
                % if you get a write permission error here, close windows explorer
                % and try again
                if frame == 1 && n_file == 1;
                    imwrite(a,savename,'tif','Compression','none','WriteMode','overwrite');%, 'Description', ImageDescriptio]);
                else
                    imwrite(a,savename,'tif','Compression','none','WriteMode','append');
                end;
                waitbar(frame/n_frames,h_waitbar,['Saving corrected file' num2str(n_file) '/' num2str(length(filename)) ' frame: ' num2str(frame) '/' num2str(n_frames)]);
            end
            close(h_waitbar);
            clear im_s_snips im_c
        end
    end
end


% --- Executes on button press in turbreg.
function turbreg_Callback(hObject, eventdata, handles)
% hObject    handle to turbreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Copied from turboreg section of HMM above + modded

%Create initial reference for parsed file (im_t)
tic
master_filelist = get(handles.filelist,'UserData');
autoRef = get(handles.AutomaticRef,'value');
parsedfile = get(handles.parsedfile,'value');
channels = get(handles.channels,'value');
multisession = get(handles.multisession,'value');

if multisession == 1;
    num_session = max(master_filelist{end}{4});
elseif multisession == 0;
    filelist = master_filelist;
    num_session = 1;
end

for session = 1:num_session
    
    if multisession == 1;
        num_session = max(master_filelist{end}{4});
    elseif multisession == 0;
        filelist = master_filelist;
        num_session = 1;
    end
    
    if multisession == 1;
        session_index = cellfun(@(x) x{4}==session,master_filelist);
        filelist = master_filelist(session_index);
    end
    
    % create parsed ref if selected
    if parsedfile == 1 && autoRef == 1;
        APMotionCorrect_parsedRef
    end
    % run turboreg
    APMotionCorrect_turboreg
    % get rid of temporary ref files
    if exist([savedir filesep 'parsedAverage.tif'],'file')
        delete([savedir filesep 'parsedAverage.tif']);
    end
    
    % create final average of all files
    
    % get image size from first file
    imageinfo=imfinfo(filelist{1}{2,1},'tiff');
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
    
    w = waitbar(0,'Creating average from corrected');
    ref_numframes = 0;
    im_avg = zeros(N,M);
    im_max = zeros(N,M);
    % get total number of frames from all files
    for pfiles2 = 1:size(filelist,1);
        ref_iminfo=imfinfo(filelist{pfiles2}{2,1},'tiff');
        ref_numframes_long=length(ref_iminfo);
        ref_numframes=ref_numframes + ref_numframes_long/channels;
    end
    for pfiles = 1:size(filelist,1);
        % Get paths/filenames from GUI (stored in filelist for each
        savedir = filelist{filenum}{1,1};
        filename = filelist{filenum}{2,1};
        reffilename = filelist{filenum}{3,1};
        [filename_path,filename_file,filename_ext] = fileparts(filename);
        % Get saved name of files
        savecorrected = [savedir filesep filename_file '_Turboreg.tif'];
        % Build up average image and max image
        ref_iminfoCurrent=imfinfo(filelist{pfiles}{2,1},'tiff');
        ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
        for ref_frame = 1:ref_numframesCurrent
            im_build_curr = double(imread(savecorrected,'tiff',ref_frame));
            im_avg(:,:) = im_avg(:,:)+im_build_curr./ref_numframes;
            im_max_temp(:,:,1) = im_build_curr;
            im_max_temp(:,:,2) = im_max;
            im_max(:,:) = max(im_max_temp,[],3);
        end
        waitbar(pfiles/size(filelist,1),w,'Creating average from corrected');
    end
    % save temporary combined average and max
    [filename_path,filename_file,filename_ext] = fileparts(filelist{filenum}{2,1});
    im_avg = uint16(im_avg);
    im_max = uint16(im_max);
    for windows_lock = 1:100
        try
            imwrite(im_avg,[savedir filesep filename_file 'AVG.tif'],'tif','Compression','none','WriteMode','overwrite');
            imwrite(im_max,[savedir filesep filename_file 'MAX.tif'],'tif','Compression','none','WriteMode','overwrite');
            break;
        catch me
            pause(0.2);
        end
    end
    clear im_t
    close(w);
    toc
end

% --- Executes on button press in parsedfile.
function parsedfile_Callback(hObject, eventdata, handles)
% hObject    handle to parsedfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parsedfile



% --- Executes on button press in multisession.
function multisession_Callback(hObject, eventdata, handles)
% hObject    handle to multisession (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multisession


% --- Executes on selection change in motionBias.
function motionBias_Callback(hObject, eventdata, handles)
% hObject    handle to motionBias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns motionBias contents as cell array
%        contents{get(hObject,'Value')} returns selected item from motionBias


% --- Executes during object creation, after setting all properties.
function motionBias_CreateFcn(hObject, eventdata, handles)
% hObject    handle to motionBias (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
