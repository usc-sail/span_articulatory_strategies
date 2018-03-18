function strategies = getStrategies(configStruct,folder,q)
% GETSTRATEGIES - get strategies
% 
% INPUT: 
%  Variable name: config_struct
%  Size: 1x1
%  Class: struct
%  Description: Fields correspond to constants and hyperparameters. 
%  Fields: 
%  - out_path: (string) path for saving MATLAB output
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

f = configStruct.f;

% Load contour data, constriction degree measurements, and factors
load(fullfile(configStruct.in_path,folder),'contour_data')

% Get weights.
X = [contour_data.X, contour_data.Y];
X = X - mean(X);
W = X*pinv(contour_data.U_gfa');
nw = size(W,2);

n = length(contour_data.tv{1}.cd); % number of video frames
nz = length(contour_data.tv);
Z = zeros(n,nz);
for i=1:nz
    % Get constriction degrees.
    Z(:,i) = contour_data.tv{i}.cd;
end

lip_idx = 1;
alv_idx = 2;
pal_idx = 3;
vel_idx = 4;
phar_idx = 5;
use_u = lip_idx*cellfun(@(x) contains(x,'_apa_'), contour_data.file_list) ...
    + alv_idx*cellfun(@(x) contains(x,'_ata_'), contour_data.file_list) ...
    + pal_idx*cellfun(@(x) contains(x,'_aia_'), contour_data.file_list) ...
    + vel_idx*cellfun(@(x) contains(x,'_aka_'), contour_data.file_list);
use = zeros(size(Z,1),1);
for i=1:length(use_u)
    use = use + use_u(i)*(contour_data.files == i);
end

% Use central difference formula to get time derivative of weights and
% constriction degrees. 
[~,dwdt] = get_grad(Z,W,1,contour_data.files);

% Get q nearest neighbors of each articulator parameter value.
fn = round(f*n);
[idx, dist] = knnsearch(W,W,'dist','euclidean','K',fn);

% Initialize containers.
jaw = zeros(n,2);
lip = zeros(n,2);
tng = zeros(n,2);
vel = zeros(n,2);

% Initialize constants.
jaw_idxs = 1:q.jaw;
tng_idxs = (q.jaw+1):(q.jaw+q.tng);
lip_idxs = (q.jaw+q.tng+1):(q.jaw+q.tng+q.lip);
vel_idxs = (q.jaw+q.tng+q.lip+1):(q.jaw+q.tng+q.lip+q.vel);
Pjaw = zeros(nw);
Pjaw(jaw_idxs,jaw_idxs) = eye(length(jaw_idxs));
Ptng = zeros(nw);
Ptng(tng_idxs,tng_idxs) = eye(length(tng_idxs));
Plip = zeros(nw);
Plip(lip_idxs,lip_idxs) = eye(length(lip_idxs));
Pvel = zeros(nw);
Pvel(vel_idxs,vel_idxs) = eye(length(vel_idxs));

for i=1:n
    if use(i) ~= pal_idx
        % Estimate the jacobian J of the forward map at point w(i,:)
        G = lscov([ones(length(idx(i,:)),1) W(idx(i,:),:)], Z(idx(i,:),use(i)), ...
            arrayfun(@(u) weight_fun(u), dist(i,:)./dist(i,end)));
        J = G(2:end)';
        
        % Determine the contributions of articulators to change in
        % constriction degree.
        jaw(i,:) = [J*Pjaw*dwdt(i,:)' NaN];
        lip(i,:) = [J*Plip*dwdt(i,:)' NaN];
        tng(i,:) = [J*Ptng*dwdt(i,:)' NaN];
        vel(i,:) = [J*Pvel*dwdt(i,:)' NaN];
    else
        % Estimate the jacobian J of the forward map at point w(i,:)
        G = lscov([ones(length(idx(i,:)),1) W(idx(i,:),:)], Z(idx(i,:),[use(i) phar_idx]), ...
            arrayfun(@(u) weight_fun(u), dist(i,:)./dist(i,end)));
        J = G(2:end,:)';
        
        % Determine the contributions of articulators to change in
        % constriction degree.
        jaw(i,:) = J*Pjaw*dwdt(i,:)';
        lip(i,:) = J*Plip*dwdt(i,:)';
        tng(i,:) = J*Ptng*dwdt(i,:)';
        vel(i,:) = J*Pvel*dwdt(i,:)';
    end
end

strategies = struct('jaw',jaw,'lip',lip,'tng',tng,'vel',vel,'tv',use,'files',contour_data.files,'file_list',contour_data.file_list);

end