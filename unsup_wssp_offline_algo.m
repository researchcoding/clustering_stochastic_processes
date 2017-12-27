function [I, dm] = unsup_wssp_offline_algo(m,m_idx,k,varargin)
% UNSURP_WSSP_OFFLINE_ALGO: perform offline clustering algorithm and with
%                           respect to observed times series, with known
%                           number of clusters
% 
% INPUT
%   m:     [matrix] observed time series matrix. Each row is an observed time series;
%   m_idx: [matrix] observed values (labeled '1') of m matrix. Size much be the
%          same as m;
%   k: [scalar] pre-defined number of clusters (kappa),
%   varargin: [matrix] use the input distance matrix
%
% OUTPUT
%   I:      group number of each observed time series labeled from 1 to k.
%           Length is same with the row number of m;
%   dm:     distance matrix (symmetric) of observed time series in m.

    if ~isequal(size(m), size(m_idx))
        error(' Size of observation matrix is not equal to size of observed index matrix');
    elseif size(m,1) < k
        error(' Number of observations is less than number of clusters.')
    end

    N = size(m, 1);
    remain_idx = 1:N;
    clustered_idx = [];
    I = zeros(N,1);
    dm = [];
    
    if k == 1 
        % trivial case
        I = ones(size(m, 1), 1);
    else 
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
                    dm(i,j) = dist_ts(m(i,com_obs_ind),m(j,com_obs_ind),'cov');
                    dm(j,i) = dm(i,j);
                end
            end
        end

        % first group
        init_idx = find(max(abs(dm)) == max(max(abs(dm))),1);
        I(init_idx) = 1;
        remain_idx(remain_idx == init_idx) = [];
        clustered_idx = [clustered_idx init_idx];
        
        % second group
        d = min(dm(init_idx, remain_idx), [], 1);
        [~,idx] = find(d == max(d),1);
        I(remain_idx(idx)) = 2;
        clustered_idx = [clustered_idx remain_idx(idx)];
        remain_idx(remain_idx == remain_idx(idx)) = [];
        
        if k >= 3
            % separate groups
            for i = 3:k
                % distance based clustering
                d = min(dm(clustered_idx(1:(i-1)), remain_idx), [], 1);
                [~,idx] = find(d == max(d),1);
                I(remain_idx(idx)) = i;
                clustered_idx = [clustered_idx remain_idx(idx)];
                remain_idx(remain_idx == remain_idx(idx)) = [];
            end
        end
        
        % assign remaining points - distance based
        for i = 1:length(remain_idx)
            d = dm(remain_idx(i), clustered_idx);
            [~,idx] = find(d == min(d), 1);
            I(remain_idx(i)) = I(clustered_idx(idx));
            clustered_idx = [clustered_idx remain_idx(i)];
        end
    end    
end 