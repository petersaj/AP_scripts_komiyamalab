%% create average image

im_s = mean(handles.im,3);
    
 %% Mukamel autosort

[fn_path fn_file fn_ext] = fileparts(tiff_fullfilename);
fn = tiff_fullfilename
outputdir = fn_path;

[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, [], nPCs, [], outputdir, []);

% if you want to manually choose PCs:
[PCuse] = CellsortChoosePCs(fn,mixedfilters);
% if you want to use all PCs:
%PCuse = 1:nPCs;

% if you want to plot PC spectrum
%CellsortPlotPCspectrum(fn,CovEvals,PCuse)

nIC = length(PCuse);
[ica_sig,ica_filters,ica_A,numiter] = CellsortICA(mixedsig,mixedfilters,CovEvals,PCuse,mu,nIC,[],termtol,maxrounds);
disp_mode = 'series';
tlims = [0 600]; %look up
ratebin = 1;
plottype = 1; % cellular signals only
ICuse = 1:nIC;
spt = [];
spc = [];
% if you want to plot all ICA traces
%CellsortICAplot(disp_mode,ica_filters,ica_sig,f0,tlims,dt,ratebin,plottype,ICuse,spt,spc)
plotting = 0; % don't show process - slows it down
[ica_segments,segmentlabel,segcentroid] = CellsortSegmentation(ica_filters,smwidth,thresh,arealims,plotting);
% if you want to plot all segmented(?) traces
%CellsortICAplot(disp_mode,ica_segments,ica_sig,f0,tlims,dt,ratebin,plottype,ICuse,spt,spc)