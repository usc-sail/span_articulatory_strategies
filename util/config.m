function configStruct = config
% CONFIG - set constants and parameters of the analysis
% 
% INPUT
%  none
% 
% FUNCTION OUTPUT:
%  Variable name: configStruct
%  Size: 1x1
%  Class: struct
%  Description: Fields correspond to constants and hyperparameters. 
%  Fields: 
%  - outPath: (string) path for saving MATLAB output
%  - aviPath: (string) path to the AVI files
%  - graphicsPath: (string) path to MATALB graphical output
%  - trackPath: (string) path to segmentation results
%  - manualAnnotationsPath: (string) path to manual annotations
%  - timestamps_file_name_<dataset>: (string) file name with path of 
%      timestamps file name for each data-set <dataset> of the analysis
%  - folders_<dataset>: (cell array) string folder names which belong to 
%      each data-set <dataset> of the analysis
%  - tasks: (cell array) string identifiers for different tasks
%  - FOV: (double) size of field of view in mm^2
%  - Npix: (double) number of pixels per row/column in the imaging plane
%  - framespersec_<dataset>: (double) frame rate of reconstructed real-time
%      magnetic resonance imaging videos in frames per second for each 
%      data-set <dataset> of the analysis
%  - ncl: (double array) entries are (i) the number of constriction 
%      locations at the hard and soft palate and (ii) the number of 
%      constriction locations at the hypopharynx (not including the 
%      nasopharynx).
%  - f: (double) hyperparameter which determines the percent of data used 
%      in locally weighted linear regression estimator of the jacobian; 
%      multiply f by 100 to obtain the percentage
%  - verbose: controls non-essential graphical and text output
% 
% SAVED OUTPUT: 
%  none
% 
% EXAMPLE USAGE: 
%  >> configStruct = config;
% 
% Tanner Sorensen
% Signal Analysis and Interpretation Laboratory
% Feb. 14, 2017

% paths
inPath = '../span_contour_processing/mat';
outPath = 'mat';
graphicsPath = 'graphics';
manualAnnotationsPath = 'manual_annotations';
trackPath = '../segmentation_results';

% timestamps file name
timestamps_file_name = 'manual_annotations/segments_ms.csv';

% array constants
folders = dir(fullfile(trackPath));
folders = {folders.name};
folders = folders(3:end);

% fixed parameters
FOV = 200; % 200 mm^2 field of view 
Npix = 68; % 68^2 total pixels
%spatRes = FOV/Npix; % spatial resolution
framespersec = 1/(2*0.006004); % compare to earlier 1/(7*6.164)*1000
ncl = [3 1];

% free parameters
f = 0.1;

% control printed output
verbose = true;

% make the struct object
configStruct = struct('inPath',inPath,'outPath',outPath,...
    'graphicsPath',graphicsPath,...
    'timestamps_file_name',timestamps_file_name,...
    'manualAnnotationsPath',manualAnnotationsPath,...
    'folders',{folders},...
    'FOV',FOV,'Npix',Npix,...%'spatRes',spatRes,...
    'framespersec',framespersec,...
    'f',f,'ncl',ncl,'verbose',verbose);

end