function [obs, cluster_ind] = sim_wssp_paths(n_clusters, obs_num_clusters, H, ...
                                             num_points_one_path, num_sim_paths)
% SIM_WSSP_PATHS: Simulate weakly stationary stochastic process (WSSP)
%                 suggested by Shields (1996). The code outlines the main
%                 simulation study steps of Peng, Rao and Zhao (2017)
% 
% INPUT
%   n_clusters: [scaler] total number of clusters;
%   obs_num_clusters: [scalar or vector] specify number of observed stochastic 
%                     processes in each cluster. If the input is a scalar,
%                     equal number of observed processes will be assumed
%                     for each cluster;
%   H: [vector] the Hurst parameter of fBm in each cluster
%   num_points_one_path: [scalar] number of observed points in each path;
%   num_sim_paths: [scalar] number of iterations in simulation study.
%   
% OUTPUT
%   obs: [matrix] the 3-dimensional matrix with 1st dimensional on each
%        individual observation, 2nd dimension on different observations,
%        and 3rd dimension on different simulation iterations;
%   cluster_ind: [vector] the cluster labels of each observation.
% 
% REFERENCE
%   

% beginning of the code
if length(obs_num_clusters) ~= 1 && obs_num_clusters ~= n_clusters
    error('The length of obs_num_clusters vector is not equal to the number of clusters. \n')
end
if length(H) ~= n_clusters
    eror('The length of alpha vector is not equal to number of clusters')
end
% initialization
if length(obs_num_clusters) == 1 && n_clusters~= 1
    obs_num_clusters = obs_num_clusters * ones(n_clusters, 1);
end
cluster_ind = ones(sum(obs_num_clusters), 1);
if n_clusters > 1
    for i = 2:n_clusters
        start_ind = sum(obs_num_clusters(1:(i-1))) + 1;
        end_ind = sum(obs_num_clusters(1:i));
        cluster_ind(start_ind:end_ind) = i * ones(obs_num_clusters(i), 1);
    end
end
obs = zeros(num_points_one_path, sum(obs_num_clusters), num_sim_paths);

% generate weakly stationary process
for i = 1:num_sim_paths
    for z = 1:sum(obs_num_clusters)
        % first order difference of fBm paths
        obs(2:end, z, i) = diff(fbmlevinson(num_points_one_path, H(cluster_ind(z))), 1);
    end
end
% end of function
end

