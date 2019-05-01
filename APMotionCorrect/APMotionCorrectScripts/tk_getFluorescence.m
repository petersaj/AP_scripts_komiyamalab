%% take motion corrected images, draw ROIs and get cell image

%% ROIs
startResponseTracker
%choose a file in the middle, draw ROIs
% (8,0.4,20,30) works pretty well
% 'Force aspect ratio to square' on and average all frames, define ROIs
% save as ----FAR.rtrk, then unclick 'Forace aspect ratio to square', save
% as ----.rtrk.
%% make result files

parentdir  = uigetdir('E:\ImagingData','Pick a date directory');
filedir = [parentdir filesep];
files = dir(filedir);
fileNames = {files.name};
Newdir = 'Result';
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(parentdir,Newdir);
h=waitbar(0,'making result files...');
for i=3:length(files)
    waitbar(i/length(files));
    [PATHSTR,NAME,EXT,VERSN] = fileparts(files(i).name);
    if isequal(EXT, '.tif')
        sourcefile=([filedir filesep files(i).name]);
        [sourcefilepath, sourcefilename] = fileparts(sourcefile);
        [PATHSTR,NAME,EXT,VERSN] = fileparts(sourcefile);
        load([sourcefilepath, NAME, '.rtrk'],'-mat')
        info=imfinfo(sourcefile);
        Result.CellImageOriginal=zeros(size(annotations,2),size(info,2));
        for i=1:size(info,1)
            im=imread(sourcefile,i);
            for j=1:size(annotations,2)
                Result.CellImageOriginal(j,i) = sum(im(annotations(j).pixels))/size(annotations(j).pixels,1);
            end
        end
        save([parentdir, filesep, Newdir, filesep, NAME, '.Result'],'Result')
    end
end
close(h)