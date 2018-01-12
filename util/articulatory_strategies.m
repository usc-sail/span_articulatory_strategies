function articulatory_strategies(configStruct)
% ARTICULATORY_STRATEGIES - get articulatory strategies, print out the 
% quantification, and make graphics
% 
% INPUT
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
%  Variable name: dataset
%  Size: arbitrary
%  Class: char
%  Description: determines which data-set to analyze; picks out the
%  appropriate constants from configStruct.
% 
% FUNCTION OUTPUT: 
%  none
% 
% SAVED OUTPUT:
%  Path: configStruct.pathOut
%  File name: strategies_<dataset>.mat
%  Variable name: strategies
%  Size: 1x1
%  Class: struct
%  Description: Struct with fields for each participant. Each field is 
%  itself a struct with the following fields.
%  - jaw: Nx6 array of double; entries are jaw contributions to change in
%  each of 6 constriction degrees (columns) in N real-time magnetic 
%  resonance imaging video frames (rows)
%  - lip: Nx6 array of double; entries are lip contributions to change in
%  each of 6 constriction degrees (columns) in N real-time magnetic 
%  resonance imaging video frames (rows)
%  - tng: Nx6 array of double; entries are tongue contributions to change 
%  in each of 6 constriction degrees (columns) in N real-time magnetic 
%  resonance imaging video frames (rows)
%  - tng: Nx6 array of double; entries are velum contributions to change 
%  in each of 6 constriction degrees (columns) in N real-time magnetic 
%  resonance imaging video frames (rows)
%  - dz: Nx6 array of double; entries are change in each of 6 constriction 
%  degrees (columns) in N real-time magnetic resonance imaging video frames
%  (rows)
%  - dw: Nx8 array of double; entries are change in 8 factor coefficients 
%  (columns) in N real-time magnetic resonance imaging video frames (rows)
%  - cl: (cell array of strings) identifier for the 6 places of 
%  articulation
% 
%  Path: configStruct.pathOut
%  File name: strategies_<dataset>.csv
%  Description: This is a spreadsheet with columns for lip, tongue, velum, 
%    and jaw biomarkers, total elapsed change in z, folder, and linguistic 
%    task 
% 
% Tanner Sorensen
% Signal Analysis and Interpretation Laboratory
% Feb. 14, 2017

outPath = configStruct.outPath;
folders = configStruct.folders;

strategies = struct;
for i=1:length(folders)
    participant =folders{i};
    disp(participant)
    strategies.(participant) = getStrategies(configStruct,participant);
    print_error(strategies.(participant))
end
save(fullfile(outPath,'strategies.mat'),'strategies');

end