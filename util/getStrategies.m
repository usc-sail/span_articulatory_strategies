function strategies = getStrategies(configStruct,folder)
% GETSTRATEGIES - get strategies
% 
% INPUT: 
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
% FUNCTION OUTPUT:
%  Variable name: strategies
%  Size: 1x1
%  Class: struct
%  Description: Struct with the following fields.
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
% SAVED OUTPUT: 
%  none
% 
% Tanner Sorensen
% Signal Analysis and Interpretation Laboratory
% Dec. 20, 2016

inPath = configStruct.inPath;
f = configStruct.f;

% Load contour data, constriction degree measurements, and factors
load(fullfile(inPath,'contourdata.mat'))
load(fullfile(inPath,'tv.mat'));
load(fullfile(inPath,'U_gfa.mat'));

% Get weights.
xy = [contourdata.(folder).X, contourdata.(folder).Y];
xy = zscore(xy);
w = xy*U_gfa.(folder);
w = zscore(w);
nw = size(w,2);

n = length(tv.(folder).tv{1}.cd); % number of video frames
nz = length(tv.(folder).tv);
z = zeros(n,nz);
use = cell(1,nz);
cl = cell(1,nz);
for i=1:nz
    % Get constriction degrees.
    z(:,i) = tv.(folder).tv{i}.cd;
    
    % The 6 entries of the cell array USE are the factor IDs which are
    % associated with the 6 constriction degree locations. 
    % factor ID key: 1-jaw, 2-5-tng, 6-7-lips, 8-velum
    cl{i} = tv.(folder).tv{i}.cl;
    if strcmp(cl{i},'bilabial')
        use{i} = [1 6 7];
    elseif strcmp(cl{i},'pharU')
        use{i} = 8;
    elseif strcmp(cl{i},'alv') || strcmp(cl{i},'pal') || strcmp(cl{i},'pharL')
        use{i} = 1:5;
    else % softpal
        use{i} = [1:5 8];
    end
end

% Use central difference formula to get time derivative of weights and
% constriction degrees. 
[dzdt,dwdt] = getGrad(z,w,1,contourdata.(folder).File);

% Get q nearest neighbors of each articulator parameter value.
q = round(f*n);
[idx, dist] = knnsearch(w,w,'dist','euclidean','K',q);

% Initialize containers.
J = zeros(nz,nw);
jaw = zeros(n,nz);
lip = zeros(n,nz);
tng = zeros(n,nz);
vel = zeros(n,nz);
dz = zeros(n,nz);
dw = zeros(n,nw);

% Initialize constants.
Pjaw = eye(nw);
Pjaw(2:end,2:end) = 0;
Ptng = eye(nw);
Ptng([1 6:end],[1 6:end]) = 0;
Plip = eye(nw);
Plip([1:5 8:end],[1:5 8:end]) = 0;
Pvel = eye(nw);
Pvel(1:7,1:7) = 0;

fprintf('['), ell=1;
for i=1:n
    dw(i,:) = dwdt(i,:)';
    for j=1:nz
        % Estimate the jacobian J of the forward map at point w(i,:)
        J(j,use{j}) = lscov(dwdt(idx(i,:),use{j}), ...
            dzdt(idx(i,:),j), ...
            arrayfun(@(u) weightfun(u), dist(i,:)./dist(i,end)));

        % Determine the contributions of articulators to change in
        % constriction degree.
        jaw(i,j) = J(j,:)*Pjaw*dw(i,:)';
        lip(i,j) = J(j,:)*Plip*dw(i,:)';
        tng(i,j) = J(j,:)*Ptng*dw(i,:)';
        vel(i,j) = J(j,:)*Pvel*dw(i,:)';
        dz(i,j) = dzdt(i,j);
    end
    
    if i/n>ell*0.05, fprintf('='), ell=ell+1; end
end
fprintf(']')

strategies = struct('jaw',jaw,'lip',lip,'tng',tng,'vel',vel,'dz',dz,'dw',dw,'cl',{cl});

end