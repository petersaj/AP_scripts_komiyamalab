%% take motion corrected images, draw ROIs and get cell image

%% ROIs

% uncomment this if response tracker isn't open 
%startResponseTracker


%choose a file in the middle, draw ROIs
% (8,0.4,20,30) works pretty well
% 'Force aspect ratio to square' on and average all frames, define ROIs
% save as ----FAR.rtrk, then unclick 'Forace aspect ratio to square', save
% as ----.rtrk.
%% make result files

[filenames parentdir] = uigetfile('*.rtrk','Pick files','MultiSelect','On');
[savedir] = uigetdir('Pick save directory');
if ~iscell(filenames);
    filenames = mat2cell(filenames);
end

h=waitbar(0,'making result files...');
for i=1:length(filenames)
    waitbar(i/length(filenames));
    sourcefile=([parentdir filenames{i}(1:end-5) '.tif']);
    load([parentdir filenames{i}],'-mat')
    info=imfinfo(sourcefile);
    Result.CellImageOriginal=zeros(size(annotations,2),size(info,2));
    for frames=1:size(info,1)
        im=imread(sourcefile,frames);
        for j=1:size(annotations,2)
            Result.CellImageOriginal(j,frames) = sum(im(annotations(j).pixels))/size(annotations(j).pixels,1);
        end
        waitbar(frames/size(info,1),h,'Creating cell trace')
    end
    save([savedir, filesep, filesep, filenames{i}(1:end-5), '.Result'],'Result')
end
close(h)