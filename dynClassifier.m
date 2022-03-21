function [subgroups_ID centroids_ID distM] = dynClassifier(X,n,varargin)
%% DYNCLASSIFIER(X,k,varargin) is designed to cluster cellular dynamics into 2n subgroups using kshape and binary tree
% INPUT:
%       X: time-series data in a matrix with each column representing a
%       time-course data
%       n: the level of binary structure, e.g n=3 means three-level binary
%       tree with total 2x2x2 = 8 subgroups
%       varargin{1}(Optional): repeat time of running clustering algorithm, 10
%       by default
%       varargin{2}(Optional): distance matrix, only when distance matrix is obtained before   
% OUTPUT:
%       subgroups_ID: matrix of cluster IDs with the xth column
%       representing the cluster ID for xth level binary tree
%       Centroids: cell of the most representative trajectories in each cluster
%       distM: distance matrix
% 
% written by Caibin Sheng(shengcaibin@gmail.com), Loewer lab, TU Darmstadt

% check input
if nargin==2
    repeat = 20;    % default repeat time is 20
elseif nargin==3
    repeat = varargin{1};
elseif nargin ==4
    repeat = varargin{1};
    distM = varargin{2}; % run only clustering algorithm if shape-based distance (SBD) matrix is already obtained
end

[~,num_of_data] = size(X); % row of X: time point; col of X: number of time-course trajectories
subgroups_ID = nan(num_of_data,n);
centroids_ID = cell(1,n);

for i =1:n
    % get the first level clusters
    if i == 1
        [cluster_ID Centroids distM] = kshape2(X,2,repeat);
        subgroups_ID(:,i) = cluster_ID;
        centroids_ID{i} = Centroids;
        
    % get clusters for 2 to even lower levels
    else
        
        subgrpsIDs = unique(subgroups_ID(:,i-1));
 
        for j = 1:length(subgrpsIDs)
            index = find(subgroups_ID(:,i-1) == subgrpsIDs(j));
            X_subgroups = X(:,index);
            
            distM_subgroups = distM(index,index);
            [cluster_ID Centroids1 ~] = kshape2(X_subgroups,2,repeat,distM_subgroups);
            subgroups_ID(index,i) = cluster_ID + 10*subgrpsIDs(j);
            centroids_ID{i}(j,:) = Centroids1;
        end
    end
end
end
