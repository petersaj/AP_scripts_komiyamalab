function varargout = APMotionCorrect_par(varargin)
% APMOTIONCORRECT MATLAB code for APMotionCorrect.fig
%      APMOTIONCORRECT, by itself, creates a new APMOTIONCORRECT or raises the existing
%      singleton*.
%
%      H = APMOTIONCORRECT_PAR_par returns the handle to a new APMOTIONCORRECT_PAR_par or the handle to
%      the existing singleton*.
%
%      APMOTIONCORRECT_PAR_par('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APMOTIONCORRECT_PAR_par.M with the given input arguments.
%
%      APMOTIONCORRECT_PAR_par('Property','Value',...) creates a new APMOTIONCORRECT_PAR_par or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before APMOTIONCORRECT_PAR_par_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to APMOTIONCORRECT_PAR_par_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help APMOTIONCORRECT_PAR_par

% Last Modified by GUIDE v2.5 17-Mar-2011 13:47:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @APMOTIONCORRECT_PAR_par_OpeningFcn, ...
    'gui_OutputFcn',  @APMOTIONCORRECT_PAR_par_OutputFcn, ...
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


% --- Executes just before APMOTIONCORRECT_PAR_par is made visible.
function APMOTIONCORRECT_PAR_par_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to APMOTIONCORRECT_PAR_par (see VARARGIN)

% Choose default command line output for APMOTIONCORRECT_PAR_par
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes APMOTIONCORRECT_PAR_par wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = APMOTIONCORRECT_PAR_par_OutputFcn(hObject, eventdata, handles)
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
tic

%Create initial reference for parsed file (im_t)
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
        session_index = cellfun(@(x) x{4}==session,master_filelist);
        filelist = master_filelist(session_index); 
    end
    
    if parsedfile == 1 && autoRef == 1;
        savedir = filelist{1}{1,1};
        
        % get image size from first file
        imageinfo=imfinfo(filelist{1}{2,1},'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        w = waitbar(0,'Turboreg: Creating Reference from Parsed');
        ref_numframes = 0;
        im_t = zeros(N,M);
        % get total number of frames from all files
        for pfiles2 = 1:size(filelist,1);
            ref_iminfo=imfinfo(filelist{pfiles2}{2,1},'tiff');
            ref_numframes_long=length(ref_iminfo);
            ref_numframes=ref_numframes + ref_numframes_long/channels;
        end
        for pfiles = 1:size(filelist,1);
            ref_iminfoCurrent=imfinfo(filelist{pfiles}{2,1},'tiff');
            ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
            for ref_frame = 1:channels:ref_numframesCurrent*channels
                im_t(:,:)=im_t(:,:)+double(imread(filelist{pfiles}{2,1},'tiff',ref_frame))./ref_numframes;
            end
            waitbar(pfiles/size(filelist,1),w,'Turboreg: Creating Reference from Parsed');
        end
        % save temporary combined average
        im_t = uint16(im_t);
        imwrite(im_t,[savedir filesep 'tempFirstAverage.tif'],'tif','Compression','none','WriteMode','overwrite');
        clear im_t
        close(w);
        
    end
    
    filelist = get(handles.filelist,'UserData');
    if exist('filelist')
        
        % loop through all files in list
        filenum = 1;
        while filenum <= size(filelist,1)
            
            % Get paths/filenames from GUI
            savedir = filelist{filenum}{1,1};
            filename = filelist{filenum}{2,1};
            reffilename = filelist{filenum}{3,1};
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
            parsedfile = get(handles.parsedfile,'value');
            
            % Do turboreg if selected
            
            % if parsed + currently autorefd: use avg reference created earlier
            if parsedfile == 1 && autoRef == 1;
                autoRef = 0;
                reffilename = [savedir filesep 'tempFirstAverage.tif'];
            end
            
            if turboreg == 1 && ~exist([savedir filesep 'tempTurboregAverage.tif'],'file');
                
                savename = [savedir filesep filename_file '_TurboregHMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
                
                imageinfo=imfinfo(filename,'tiff');
                numframes=length(imageinfo);
                M=imageinfo(1).Width;
                N=imageinfo(1).Height;
                %im_s=zeros(N,M,numframes/channels);
                
                % Autoref: make avg image, temp zeros for source image
                if autoRef == 1;
                    w = waitbar(0,'Creating average image');
                    im_t = zeros(N,M);
                    for u = 1:channels:numframes
                        im_temp(:,:)=double(imread(filename,'tiff',u));
                        im_t(:,:) = im_t + im_temp./(numframes/channels);
                        waitbar(u/numframes,w);
                    end
                    im_s = zeros(N,M);
                end
                
                % Get reference or create average for reference
                if autoRef == 0;
                    im_t = imread(reffilename,'tiff');
                end;
                
                % If you've got NaN, make 0 (keep in mind for correlation)
                %im_s(im_s<1)=0;
                
                % Add ImageJ java controls for turboreg
                spath = javaclasspath('-static');
                spath = cell2mat(spath(1));
                javafolder = strfind(spath,['java' filesep]);
                javafolder = spath(1:javafolder+3);
                
                javaaddpath([javafolder filesep 'jar' filesep 'ij.jar'])
                javaaddpath([javafolder filesep 'jar' filesep 'mij.jar'])
                javaaddpath([javafolder])
                
                if autoRef == 0;
                    Turboreg_AP(reffilename,filename,savename,channels);%(targetfilename, sourcefilename, savefilename, channels)
                elseif autoRef == 1;
                    Turboreg_AP(im_t,filename,savename,channels);
                end
                clear im_s
                clear im_t
            end
            
            
            % Use the corrected file for HMM if autoRef selected
            filename = [savedir filesep filename_file '_TurboregHMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
            % get real autoRef value again
            autoRef = get(handles.AutomaticRef,'value');
            
            if autoRef == 1 && parsedfile == 0;
                imageinfo=imfinfo(filename,'tiff');
                numframes=length(imageinfo);
                M=imageinfo(1).Width;
                N=imageinfo(1).Height;
                im_t = zeros(N,M);
                
                w = waitbar(0,'Creating reference from turboreg corrected images');
                
                for u = 1:channels:numframes
                    im_temp(:,:)=imread(filename,'tiff',u);
                    im_t = im_t + im_temp./(numframes/channels);
                    waitbar(u/numframes,w);
                end
                close(w);
                reffilename = im_t;
                % if parsedfile/autoref, do all turboreg first then take avg
            elseif autoRef == 1 && parsedfile == 1;
                if filenum ~= size(filelist,1) && ~exist([savedir filesep 'tempTurboregAverage.tif'],'file')
                    filenum = filenum+1;
                    continue
                elseif filenum == size(filelist,1) && ~exist([savedir filesep 'tempTurboregAverage.tif'],'file')
                    w = waitbar(0,'HMM: Creating Reference from Parsed Turboreg');
                    ref_numframes = 0;
                    im_t = zeros(N,M);
                    for rfiles2 = 1:size(filelist,1);
                        [tempP,tempTurboFilename_file,tempE] = fileparts(filelist{rfiles2}{2,1});
                        tempTurboFilename = [savedir filesep tempTurboFilename_file '_TurboregHMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
                        ref_iminfo=imfinfo(tempTurboFilename,'tiff');
                        ref_numframes_long=length(ref_iminfo);
                        ref_numframes=ref_numframes + ref_numframes_long;
                    end
                    for rfiles = 1:size(filelist,1);
                        [tempP,tempTurboFilename_file,tempE] = fileparts(filelist{rfiles}{2,1});
                        tempTurboFilename = [savedir filesep tempTurboFilename_file '_TurboregHMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
                        ref_iminfoCurrent=imfinfo(tempTurboFilename,'tiff');
                        ref_numframesCurrent=length(ref_iminfoCurrent);
                        for ref_frame = 1:ref_numframesCurrent
                            im_t(:,:)=im_t(:,:)+double(imread(filelist{rfiles}{2,1},'tiff',ref_frame))./ref_numframes;
                        end
                        waitbar(rfiles/size(filelist,1),w,'HMM: Creating Reference from Parsed Turboreg');
                    end
                    % save temporary combined average
                    im_t = uint16(im_t);
                    imwrite(im_t,[savedir filesep 'tempTurboregAverage.tif'],'tif','Compression','none','WriteMode','overwrite');
                    filenum = 1;
                    continue
                elseif filenum <= size(filelist,1) && exist([savedir filesep 'tempTurboregAverage.tif'],'file')
                    autoRef = 0;
                    reffilename = [savedir filesep 'tempTurboregAverage.tif'];
                end
            end
            
            % make filename for future saved HMM files
            if turboreg == 1;
                savecorrected = [savedir filesep filename_file '_TurboregHMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
            elseif turboreg == 0;
                savecorrected = [savedir filesep filename_file '_HMM_' '_MaxDx_' num2str(maxdx) '_Lambda_' get(handles.Lambda,'string')];
            end
            
            % Process with chosen reference picture, or use average
            % If turboreg was done with autoref, use new ref created above
            
            % turboreg corrected already only uses green channel
            if turboreg == 1;
                channels = 1;
            end
            
            split_createPI_ManRef(filename,savecorrected,reffilename,0,maxdx,maxdy,channels,lambda);
            load([savecorrected '_PI.mat']);
            
            
            % Calculate lambda or use given lambda
            if autoLambda == 0;
                lambda = str2num(get(handles.Lambda,'string'));
                markov_on_PIsave
            elseif autoLambda == 1;
                maxlambda;
            end
            
            % Create motion corrected image
            [fixeddata,countdata]=tk_playback_markov2(filename,savecorrected,offsets,edgebuffer,playback,lambda);
            
            % I'm not sure why I put this in here, the last function already
            % does this. I'll keep it in case it was useful for seomthing.
            
            %         timeHMM = toc;
            %         im_c = permute(fixeddata,[2 3 1]);
            %         h_waitbar = waitbar(0,'Saving corrected file');
            %         for frame = 1:length(im_c,3);
            %             a=uint16(round(im_c(:,:,frame)));
            %             % if you get a write permission error here, close windows explorer
            %             % and try again
            %             if frame == 1 && n_file == 1;
            %                 imwrite(a,savecorrected,'tif','Compression','none','WriteMode','overwrite');%, 'Description', ImageDescriptio]);
            %             else
            %                 imwrite(a,savecorrected,'tif','Compression','none','WriteMode','append');
            %             end;
            %             waitbar(frame/n_frames,h_waitbar,['Saving corrected file' num2str(n_file) '/' num2str(length(filename)) ' frame: ' num2str(frame) '/' num2str(n_frames)]);
            %         end
            
            % Delete extra files, rename finalized one
            delete([savecorrected '.tif'],[savecorrected '_PI.mat']);
            if exist([savecorrected '.mat'],'file')
                delete([savecorrected '.mat']);
            end
            if exist([savedir filesep filename_file '.mat'],'file');
                delete([savedir filesep filename_file '.mat'])
            end
            movefile([savecorrected 'FINAL.tif'],[savecorrected '.tif'],'f');
            disp(['Finished correcting: ' savecorrected])
            filenum = filenum+1;
        end
    end
    
    if exist([savedir filesep 'tempFirstAverage.tif'],'file')
        delete([savedir filesep 'tempFirstAverage.tif'], [savedir filesep 'tempTurboregAverage.tif']);
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
tic
% hObject    handle to turbreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Copied from turboreg section of HMM above + modded

%Create initial reference for parsed file (im_t)
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
        session_index = cellfun(@(x) x{4}==session,master_filelist);
        filelist = master_filelist(session_index); 
    end
    
    if parsedfile == 1 && autoRef == 1;
        savedir = filelist{1}{1,1};
        
        % get image size from first file
        imageinfo=imfinfo(filelist{1}{2,1},'tiff');
        M=imageinfo(1).Width;
        N=imageinfo(1).Height;
        
        w = waitbar(0,'Turboreg: Creating Reference from Parsed');
        ref_numframes = 0;
        im_t = zeros(N,M);
        % get total number of frames from all files
        for pfiles2 = 1:size(filelist,1);
            ref_iminfo=imfinfo(filelist{pfiles2}{2,1},'tiff');
            ref_numframes_long=length(ref_iminfo);
            ref_numframes=ref_numframes + ref_numframes_long/channels;
        end
        for pfiles = 1:size(filelist,1);
            ref_iminfoCurrent=imfinfo(filelist{pfiles}{2,1},'tiff');
            ref_numframesCurrent=length(ref_iminfoCurrent)/channels;
            for ref_frame = 1:channels:ref_numframesCurrent*channels
                im_t(:,:)=im_t(:,:)+double(imread(filelist{pfiles}{2,1},'tiff',ref_frame))./ref_numframes;
            end
            waitbar(pfiles/size(filelist,1),w,'Turboreg: Creating Reference from Parsed');
        end
        % save temporary combined average
        im_t = uint16(im_t);
        imwrite(im_t,[savedir filesep 'tempFirstAverage.tif'],'tif','Compression','none','WriteMode','overwrite');
        clear im_t
        close(w);
        
    end
    
    filelist = get(handles.filelist,'UserData');
    if exist('filelist')
        
        % loop through all files in list
        filenum = 1;
        matlabpool
        parfor filenum = 1:size(filelist,1)
            
            % Get paths/filenames from GUI
            
            %savedir = filelist{filenum}{1,1};
            filename = filelist{filenum}{2,1};
            reffilename = filelist{filenum}{3,1};
            [filename_path,filename_file,filename_ext] = fileparts(filename);
            
            %filename = get(handles.TIFFFile,'string');
            savedir = get(handles.savedir,'string');
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
            parsedfile = get(handles.parsedfile,'value');
            
            % Do turboreg if selected
            
            % if parsed + currently autorefd: use avg reference created earlier
            if parsedfile == 1 && autoRef == 1;
                autoRef = 0;
                reffilename = [savedir filesep 'tempFirstAverage.tif'];
            end
            
            if ~exist([savedir filesep 'tempTurboregAverage.tif'],'file');
                
                savename = [savedir filesep filename_file '_Turboreg'];
                
                imageinfo=imfinfo(filename,'tiff');
                numframes=length(imageinfo);
                M=imageinfo(1).Width;
                N=imageinfo(1).Height;
                %im_s=zeros(N,M,numframes/channels);
                
                % Autoref: make avg image, temp zeros for source image
                if autoRef == 1;
                    w = waitbar(0,'Creating average image');
                    im_t = zeros(N,M);
                    for u = 1:channels:numframes
                        im_temp(:,:,filenum)=double(imread(filename,'tiff',u));
                        im_t(:,:) = im_t + im_temp(:,:,filenum)./(numframes/channels);
                        waitbar(u/numframes,w);
                    end
                    im_s = zeros(N,M);
                end
                close(w);
                
                % Get reference or create average for reference
                if autoRef == 0;
                    im_t = imread(reffilename,'tiff');
                end;
                
                % If you've got NaN, make 0 (keep in mind for correlation)
                %im_s(im_s<1)=0;
                
                % Add ImageJ java controls for turboreg
                spath = javaclasspath('-static');
                spath = cell2mat(spath(1));
                javafolder = strfind(spath,['java' filesep]);
                javafolder = spath(1:javafolder+3);
                
                javaaddpath([javafolder filesep 'jar' filesep 'ij.jar'])
                javaaddpath([javafolder filesep 'jar' filesep 'mij.jar'])
                javaaddpath([javafolder])
                
                if autoRef == 0;
                    Turboreg_AP(reffilename,filename,savename,channels);%(targetfilename, sourcefilename, savefilename, channels)
                elseif autoRef == 1;
                    Turboreg_AP(im_t,filename,savename,channels);
                end
                clear im_s
                clear im_t
            end
            
            % Delete extra files, rename finalized one
            disp(['Finished correcting: ' savename])
            %filenum = filenum+1;
        end
    end
    
    if exist([savedir filesep 'tempFirstAverage.tif'],'file')
        delete([savedir filesep 'tempFirstAverage.tif'], [savedir filesep 'tempTurboregAverage.tif']);
    end
end
matlabpool close
toc

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
