% A quick overview of ICP 
% Run this for paths
addpath('ecc');
addpath('FastICA_25');
cd('GUILayout-v1p13');
install();
cd('..');

%% Data files
% set the path to your data
file = '/usr/local/lab/People/Andy/Data/AP107/140807/summed_movie/140807_AP107_summed_50.tif';

%% Load a data
im = tifread(file);
%% ImageComponentParser
% Open an icp
icp = ImageComponentParser();
% Set the raw data
icp.set_data(im);

% You can also do the same thing like this
% icp = ImageComponentParser([], im);

%   So the idea is that the analysis is broken into stages: Pre-processing,
% PCA, then ICA. And each stage you can look at various things.
% There are tabs at the bottom that control the display and give you access
% to the run buttons. You can just click run ica and every stage will get
% run sequentially if the previous stages haven't been run.
%   When you have the PCA or ICA tabs open the viewer shows you the weights
% of each components. You can click through with the slider at the bottom
% and view different components. The Component plot will show you the
% current component.
%   There is also a 3D plot so you can look at the components in that
% space. You can change which components you want to look at, and theres a
% rotation button (for help visualizing). The 3D plot will change depending
% if you have the PCA or ICA tab open.
%   Once the ICA analyis is done you can click on calc rois button in the 
% segment tab and it will come up with rois for each component (see
% calc_rois_from_components). You can also make your own ROIs and different
% ROI sets. Just use the mouse on the image to manipulate, and create ROIs.
% press ctrl+delete to delete the currently selected rois. You can add or
% delete roi sets with the + and - buttons below the roi set list. The Rois
% are defined as ovals - [x_center, y_center, x_width, y_width, angle].
% There is some infrastructure in the RoiEditor to use pixel masks, but I
% only use ovals, so this isn't that developed.
%   This is still very much a work in progress, so forgive me for bugs and
% things that are incomplete.

%% ICA Analysis

% So, you can use the gui to do this, or you can do it programatically

% First change some settings, this is what we used on your data
icp.settings.preprocessing.smooth_window = [10,10,1]; % smoothing in [m, n, t] dimensions
icp.settings.preprocessing.down_sample = [4, 4, 1];

icp.settings.ica.which_pcs = 1:200; % This means we'll use the first 200 PCs for the ICA analysis

%%
% ICP can run several different kinds of ICA, but the only thing inclued
% in the package is fastica. If you want to use 'CellsortICA' which is the
% program from that Schnitzer paper, then put it on the path and you can
% get icp to call it with this:
%icp.settings.ica.ica_func='CellsortICA';
% and there is the mu value in the settings:
%icp.settings.ica.mu=0.2;
% Check out run_ica for more info. I was trying several different kinds of
% ICA programs, but they didn't seem much better and typically took far too
% long. I've been happy with fastica

% To set the initial guess for A (the mixing matrix), there is a variable
% in the settings, icp.settings.ica.init_guess.
% init_guess == -1 means use the previously computed mixing matrix if possible
% init_guess == 0 means just use the default random
% init_guess == matrix means use the matrix as the init guess, it should yell at you if it is incorrectly formatted.

% Run each stage:
icp.run_preprocessing();
icp.run_pca();
icp.run_ica();

% Calcualte the rois from the ic components
rois = icp.calc_rois();

%% Spatial guess
% To set a spatial guess, you pass a series of masks into
% set_ica_spatial_guess. So say you ran ICA from one trial, and want to get
% components from similar location in another trial. Start by taking the
% im_data from the first trial:
im_data = icp.ica.im_data;
% Then load another trial or open another icp:
icp2 = ImageComponentParser([], tifread('the_other_filename.tif'));
% You would then set all of the settings that you would want.

% And then you should be able to just call:
icp2.set_spatial_guess(im_data);
% Which will create the appropriate init_guess matrix, and then you just
% run ICA with that init_guess
icp2.run_ica();
% Works fairly well, but not perfect.


% Thats pretty much it. Hopefully the code is sufficiently documented, but
% just let me know if you have any questions. You may want to look at
% RoiEditor, AxesEventNotifier, FigureEventNotifier for more ideas about
% how the gui works. FigureEventNotifier and AxesEventNotifier use the
% UserData field of the figures and axes, so could interfere with things if
% you use that or GUIDE.


