function [obs, cluster_ind] = sim_wssp_paths_nonGaus(n_clusters, obs_num_clusters, a, ...
                                             num_points_one_path, num_sim_paths)
% SIM_WSSP_PATHS_NONGAUS: Simulate non-Gauss weakly stationary stochastic process (WSSP)
%                 suggested by Sen and Singer (1993). The code outlines the main
%                 simulation study steps of Peng, Rao and Zhao (2018)
% 
% INPUT
%   n_clusters: [scaler] total number of clusters;
%   obs_num_clusters: [scalar or vector] specify number of observed stochastic 
%                     processes in each cluster. If the input is a scalar,
%                     equal number of observed processes will be assumed
%                     for each cluster;
%   a: [vector] the coeffcient of process X(t) = a * X(t-1) + Z(t).
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
if length(a) ~= n_clusters
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
        % simulate the process X(t) = a * X(t-1) + Z(t);
        % Z(t) = sqrt(2) * cos(tU), where U is random uniform.
        Zt = sqrt(2) * cos(linspace(0, 1, num_points_one_path).* ...
            unifrnd(0, 1, [1, num_points_one_path]));
        for j = 3:num_points_one_path
            obs(j, z, i) = a(cluster_ind(z)) * obs(j - 1, z, i) - a(cluster_ind(z))^2 * obs(j - 2, z, i) + Zt(j);
        end
    end
end
% end of function
end

