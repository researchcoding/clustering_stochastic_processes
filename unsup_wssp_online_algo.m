function [I,dm]=unsup_wssp_online_algo(m,m_idx,k,varargin)
% UNSURP_WSSP_ONLINE_ALGO: perform online clustering algorithm and with
%                          respect to observed times series, with known
%                          number of clusters
% 
% INPUT: 
%   m:      observation time series. Each row is an observed time series;
%   m_idx:  observed values (labeled '1') of m matrix. Size much be the
%           same as m;
%   k:      pre-defined number of clusters (kappa);
%   varargin: [matrix] use the input distance matrix
%
% OUTPUT:
%   I:      group number of each observed time series labeled from 1 to k.
%           Length is same with the row number of m;
%   dm:     distance matrix (symmetric) of observed time series in m.

    if ~isequal(size(m), size(m_idx))
        error(' Size of observation matrix is not equal to size of observed index matrix');
    elseif size(m,1) < k
        error(' Number of observations is less than number of clusters.')
    end
    
    N = size(m, 1);
    I = zeros(N,1);
    dm = [];
    
    if k == 1 
        % trivial case
        I = ones(N, 1);
    else
        eta = 0;
        w = zeros(1, N-k+1);
        gamma = zeros(1, N-k+1);
        mu = zeros(N-k+1,k);
        
        if nargin > 3
            dm = varargin{1};
            if ~ismatrix(dm) || size(dm, 1) ~= size(m,1), error('Invalid distance matrix'); end
        else
            % full distance matrix
            dm = zeros(N,N);
            for i = 2:N
                for j = 1:(i-1)
                    com_obs_ind = m_idx(i,:) & m_idx(j,:);
                    % covariance matrix based distance measure
                    dm(i,j) = dist_ts_log(m(i,com_obs_ind),m(j,com_obs_ind),'cov');
                    dm(j,i) = dm(i,j);
                end
            end
        end
        
        observed_points = sum(m_idx, 2);
        [~, sorted_idx] = sort(observed_points, 'descend');

        % generate N(t) - k + 1 candidate cluster centers
        for j = k:N
            selected_idx = sorted_idx(1:j);
            selected_obs = m(selected_idx, :);
            selected_obs_ind = m_idx(selected_idx, :);
            I_selected = ...
                unsup_wssp_offline_algo(selected_obs, selected_obs_ind, k, dm(selected_idx,selected_idx));
            y = [];
            for z = 1:k
                idx = selected_idx(I_selected == z);
                [~, ind] = max(sum(m_idx(idx, :),2));
                mu(j-k+1,z) = idx(ind);
                for z2 = 1:k
                    if z2 ~= z
                        idx2 = selected_idx(I_selected == z2);
                        [~, ind2] = max(sum(m_idx(idx2, :),2));
                        y = [y dm(idx2(ind2), mu(j-k+1,z))];
                    end
                end
            end
            gamma(z-k+1) = min(y);
            w(z-k+1) = 1 / (j*(j+1));
            eta = eta + gamma(z-k+1) * w(z-k+1);
        end
        
        % assign points to clusters
        for i = 1:N
            clu = zeros(1,k);
            for c = 1:k
                for j = k:N
                    clu(c) = clu(c) + gamma(z-k+1) * w(z-k+1) * dm(i, mu(j-k+1,c));
                end
            end
            [~, I(i)] = min(clu);
        end
    end    
end 