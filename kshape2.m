function [cluster_ID,Centroids,distM] = kshape2(X,k,varargin)
%% KSHAPE2(X,k,varargin) is designed to cluster time-series data X into k groups based on the shapes
% INPUT:
%       X: time-series data in a matrix with each column representing a
%       time-course data
%       k: the number of clusters
%       varargin{1}(Optional): repeat time of running clustering algorithm, 10
%       by default
%       varargin{2}(Optional): distance matrix, only when distance matrix is obtained before   
% OUTPUT:
%       cluster_ID: cluster ID
%       Centroids: the most representative time-course data in each cluster
%       distM: distance matrix 
%   
% reference:  Paparrizos, J., & Gravano, L. (n.d.). k-Shape : Efficient and Accurate Clustering of Time Series.
% written by Caibin Sheng(shengcaibin@gmail.com), Loewer lab, MDC Berlin / TU Darmstadt

% check input
if nargin==2
    repeat = 20; % default repeat time is 20
elseif nargin==3
    repeat = varargin{1};
elseif nargin ==4
    repeat = varargin{1};
    SBDmat = varargin{2}; % run only clustering algorithm if shape-based distance (SBD) matrix is already obtained
end

[time_points,num_of_data] = size(X); % row of X: time point; col of X: number of time-course trajectories


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute distance matrix
if nargin<4;
    SBDmat = NaN(num_of_data,num_of_data);


    normalized_data = zscore(X);    % z-normalization to remove scaling variance

% divide the data into small groups to compute the distance matrix in order to
% avoid overloading the PC computing power 
    size_of_subgrp = 1000; % maximum of subgroup size
    num_submat = ceil(num_of_data/size_of_subgrp);

% divide the data
    normalized_data_subgrps = mat2cell(normalized_data,[time_points],[size_of_subgrp*ones(1,num_submat-1),rem(num_of_data,size_of_subgrp)]);
% calculate shape-based distance matrix between any two divided subgroups
    SBD_subgrps = cell(num_submat,num_submat); 
    for i = 1:num_submat
        tic
        for j = 1:num_submat
            SBD_subgrps{i,j}=matxcorr(normalized_data_subgrps{1,i},normalized_data_subgrps{1,j});          
        end
        toc
        display(strcat(['-- Estimated remaining min: ',1,num2str(toc*(num_submat-i)/60),1]));
    end
    SBDmat = cell2mat(SBD_subgrps);
end
distM = SBDmat;             % return distance matrix
display ('Dissimilarity matrix done.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run clustering
display ('...........................................')
display ('Starting clustering algorithm...')

IDX_reps = cell(1,repeat);
cost_function = zeros(1,repeat);
for j = 1:repeat
    cluster_seeds = seeds(SBDmat,k); % careful seeding algorithm, see function below
    
    iteration = 0;
    IDX_reps{j} = zeros(num_of_data,1);
    IDX = ones(num_of_data,1);

    while ~all(IDX_reps{j} == IDX) && iteration<0.25*num_of_data    
        IDX_reps{j} = IDX;
        [mindist, IDX] = min(cluster_seeds,[],2); % assign each trajectory to clusters
        cost_function(j) = sum(mindist.*mindist); % get cost function
    
    % Update centroids
        for i = 1:k
            SBD_grp = SBDmat(:,IDX==i);
            if isempty(SBD_grp)
                error('Error. \n cluster number %d may be too big: reduce it or run again.',k)
            end
            ind_cen_grp = shapeExtraction (SBD_grp(IDX==i,:)); % find the relative index of the new centroids in each clusters
            cluster_seeds(:,i) = SBD_grp(:,ind_cen_grp(1)); % find and update the new controids for each clusters
        end
        iteration = iteration + 1;
        display(['    The ',num2str(j),'th experiment:',' iterations (',num2str(iteration),') energy (',num2str(cost_function(j)),')']),
    end
    % Find the index of the most representative trajectories for each clusters
    for i = 1:k
        centroids_grp = find(all(bsxfun(@eq,SBDmat,cluster_seeds(:,i))));
        loc{j}(i) = centroids_grp(1);

    end
    display (['--The sum of squared distances: ',num2str(cost_function(j)),'  The centroids: ', num2str(loc{j})])
    % Get the most representative trajectories for each clusters
    
end
    [a b] = min(cost_function);
    Centroids = loc{b};
    cluster_ID = IDX_reps{b};     
    display(['*************************************************************************************'])
    display(['After running ', num2str(repeat),' times:'])
    display (['   The sum of squared distances converges to ',num2str(a)])
    display (['   The most representative trajectories are: ', num2str(Centroids)])
end

%--------------------------------------------------------------------------
% shape extraction: find the center (the most representative trajectory) in subgroups
function IDX_centroid = shapeExtraction(data_in_subgrps)
    [~, IDX_centroid] = min(dot(data_in_subgrps,data_in_subgrps,1));
end

%--------------------------------------------------------------------------
function dist_2_subgrps = matxcorr(A,B)
% get cross-correlation between two divided subgroups (A & B) of data X
    [~,a2] = size(A);
    [~,b2] = size(B);
    [J,I] = meshgrid(1:b2,1:a2);
    dist_2_subgrps = arrayfun(@(i,j)xco(A,B,i,j),I,J);
end

%--------------------------------------------------------------------------
function dist_a_pair = xco(A,B,i,j)
% get the dissimilarity measure of pairwise trajectories {A(:,i),B(:,j)} using cross-correlation function 
    a = size(A,1);
    [c,~] = myxcorr(A(:,i),B(:,j),floor(a/2),'coeff');
    dist_a_pair = 1 - max(c);
end