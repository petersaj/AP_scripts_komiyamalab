[fn_file fn_path] = uigetfile('*.tif')
fn = [fn_path fn_file]
outputdir = uigetdir
nPCs = 90; %number of PCs

[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, [], nPCs, [], outputdir, []);
[PCuse] = CellsortChoosePCs(fn,mixedfilters);
CellsortPlotPCspectrum(fn,CovEvals,PCuse)
mu = 0.1;
nIC = length(PCuse);
termtol = 0.1;
maxrounds = 10;
[ica_sig,ica_filters,ica_A,numiter] = CellsortICA(mixedsig,mixedfilters,CovEvals,PCuse,mu,nIC,[],termtol,maxrounds);
disp_mode = 'series';
[avgfile avgpath] = uigetfile('*.tif','Choose avg file');
f0 = imread([avgpath avgfile],'tiff');
tlims = [0 600];
dt = 1;
ratebin = 3;
plottype = 1; % cellular signals only
ICuse = 1:nIC;
spt = [];
spc = [];
CellsortICAplot(disp_mode,ica_filters,ica_sig,f0,tlims,dt,ratebin,plottype,ICuse,spt,spc)

smwidth = 1.5; % std of gaussian smoother (px)
thresh = 2; % in std
arealims = [20 1000]; % min/max area (px)
plotting = 1; % show filters
[ica_segments,segmentlabel,segcentroid] = CellsortSegmentation(ica_filters,smwidth,thresh,arealims,plotting);


CellsortICAplot(disp_mode,ica_segments,ica_sig,f0,tlims,dt,ratebin,plottype,ICuse,spt,spc)
