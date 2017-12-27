function r = misclassify_rate(idx_algo, idx_true)
    N = length(idx_true);
    cluster = unique(idx_true);
    n_cluster = length(cluster);
    perm_cluster = perms(1:n_cluster);
    r_perm = zeros(size(perm_cluster,1), 1);
    for i = 1:length(r_perm)
        idx_temp = idx_true;
        for k = 1:n_cluster
            idx_temp(idx_true == k) = perm_cluster(i,k);
        end
        r_perm(i) = sum((idx_algo-idx_temp)~=0)/N;
    end
    r = min(r_perm);
end

