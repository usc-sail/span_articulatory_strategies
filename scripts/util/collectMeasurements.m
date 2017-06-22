function [jawSum, lipSum, tngSum, velSum, dzSum] = collectMeasurements( configStruct, folder, strategies )
%COLLECTMEASUREMENTS - collect articulatory strategy biomarkers for the
%right articulators in the right linguistic tasks
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
%  Variable name: folder
%  Size: arbitrary
%  Class: char
%  Description: determines which participant/scan of the data-set to 
%    analyze
% 
%  Variable name: colName
%  Size: arbitrary
%  Class: char
%  Description: determines the linguistic task (e.g., 'aia-pal-clo') which
%    is to be analyzed
%  
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

% initialize
verbose = configStruct.verbose;
strategies = strategies.(folder);
outPath = configStruct.outPath;
graphicsPath = configStruct.graphicsPath;

jaw = strategies.jaw;
lip = strategies.lip;
tng = strategies.tng;
vel = strategies.vel;
dz = strategies.dz;

% Load timestamps file
timestamps_file_name = configStruct.timestamps_file_name;
inPath = configStruct.inPath;
load(fullfile(inPath,'contourdata.mat'))

% Load constriction location names corresponding to each constriction degree.
cl_list = cell(1,length(strategies.cl));
for i=1:length(cl_list), cl_list{i} = strategies.cl{i}; end

% Load the time-stamps and the associated filename and speaker IDs. 
tab = readtable(timestamps_file_name,'Delimiter',',');
taskVar = table2cell(tab(:,5));
nFile = size(tab,1);

% Collect measurements by item.
jawSum = zeros(1,nFile);
lipSum = zeros(1,nFile);
tngSum = zeros(1,nFile);
velSum = zeros(1,nFile);
dzSum = zeros(1,nFile);

% Open file to print results
fOut = fopen(fullfile(outPath,'strategies.csv'),'w');
fprintf(fOut,'folder,task,lip,tng,jaw,vel,z\n');

% Split the files into unique file sets (indicated by a unique identifier
% in column 6 of the timestamps file)
fileSets = unique(table2cell(tab(:,6)));
for i=1:length(fileSets)
    % Determine which files from contourdata belong to the current file set
    fileMatch = find(cellfun(@(x) ~isempty(strfind(x,fileSets{i})),contourdata.(folder).fl));
    idx1 = ismember(contourdata.(folder).File,fileMatch);
    
    % Determine which constriction location to compute articulator
    % contributions towards
    idx2 = cellfun(@(x) strcmp(x,taskVar{i}),cl_list);
    
    % Compute articulator contributions
    jawSum(i) = sum(jaw(idx1,idx2));
    lipSum(i) = sum(lip(idx1,idx2));
    tngSum(i) = sum(tng(idx1,idx2));
    velSum(i) = sum(vel(idx1,idx2));
    
    % Compute total constriction degree change
    dzSum(i) = sum(dz(idx1,idx2));
    
    % Print the articulator contributions and total constriction degree change to file
    fprintf(fOut,'%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f\n',...
        folder,taskVar{i},lipSum(i),tngSum(i),jawSum(i),velSum(i),dzSum(i));
    
    % Plot the contributions of jaw and lips for up to 10 file sets
    if verbose
        figure
        bar(cumsum([jaw(idx1,idx2) tng(idx1,idx2) lip(idx1,idx2)],1),'stacked')
        nPt = sum(idx1);
        xlim([0 nPt+1]), set(gca,'XTick',[])
        pbaspect([2 1 1])
        legend({'jaw','tongue','lips'})
        print(fullfile(graphicsPath,sprintf('%s_%s_%s.png',folder,taskVar{i},fileSets{i})), '-dpng')
    end
end
fclose(fOut);

end