function clustering_kmeans(delay_data_exp, delay_data_ctrl, resdir)
%CLUSTERING_KMEANS  Perform PCA and K-means clustering on delay-period neural activity.
%   CLUSTERING_KMEANS(DELAY_DATA_EXP, DELAY_DATA_CTRL, RESDIR) performs principal component
%   analysis followed by K-means clustering on neural activity data from experimental and
%   control conditions during the delay period. The function reduces each dataset to three
%   principal components, applies K-means clustering with three clusters to the experimental
%   data using 50 replicates, and then assigns control data points to the centroids obtained
%   from the experimental data without re-clustering the control set. Cluster labels, centroids,
%   and PCA scores are saved to a MAT-file named 'clustering_res.mat' in the specified results
%   directory. A 3D scatter plot of the experimental PCA scores colored by cluster assignment,
%   with cluster centroids overlaid, is generated and saved as 'kmeanssplot.jpg'.
%
%   DELAY_DATA_EXP is an N-by-M matrix containing delay-period neural activity for N
%   experimental neurons across M features or time bins. DELAY_DATA_CTRL is a P-by-M matrix
%   for P control neurons, and must have the same number of columns M as DELAY_DATA_EXP.
%   RESDIR is a character vector specifying the output directory, which is created if it
%   does not already exist.
%
%   The saved MAT-file 'clustering_res.mat' contains the following variables:
%   clsE, the cluster labels for experimental neurons (N-by-1); clsC, the cluster labels
%   for control neurons (P-by-1) assigned to the experimental centroids; CE, the 3-by-3
%   matrix of cluster centroids from the experimental data; scoreE, the N-by-3 PCA scores
%   for the experimental data; and scoreC, the P-by-3 PCA scores for the control data.
%
%   Example:
%       clustering_kmeans(delay_act_exp, delay_act_ctrl, 'results/clustering/');
%
%   Note that the control dataset is not independently clustered. Each control neuron is
%   assigned to the nearest centroid derived from the experimental dataset.
%
%   See also PCA, KMEANS, PDIST2, SCATTER3.

%   Malek Aouadi
%   Laboratory of Systems Neuroscience, Institute of Experimental Medicine
%   Budapest, Hungary
%   2025

if ~isfolder(resdir)
    mkdir(resdir)
end

% PCA 
n_components = 3;
[~,scoreE,~] = pca(delay_data_exp, 'NumComponents', n_components);
[~,scoreC,~] = pca(delay_data_ctrl, 'NumComponents', n_components);

% Clustering
K = 3; % number of clusters

%  Cluster Dataset 1
[clsE, CE] = kmeans(scoreE, K, 'Replicates', 50);  % C1 = K x D centroids

%Assign Dataset 2 points to Dataset 1 s centroids (no re-clustering)
% Compute Euclidean distances from each point in dataset2 to each centroid in C1
distances = pdist2(scoreC, CE);  % N2 x K matrix

% Assign each point to closest centroid
[~, clsC] = min(distances, [], 2);  % labels2 is N2 x 1
save([resdir '\clustering_res.mat'],'clsE','CE', 'clsC','scoreC','scoreE','-mat');

% Scatter plot of data points colored by cluster assignment
figure;
scatter3(score(:,1), score(:, 2),score(:,3), 36, idx, 'filled');
title('K-means Clustering Results');
hold on;
scatter(C(:, 1), C(:, 2), 100, 'k', 'filled');
legend('Data Points', 'Cluster Centroids');
hold off
f='kmeanssplot';
saveas(gcf, f, 'jpg');
