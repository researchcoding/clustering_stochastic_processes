function d = dist_ts_log(A, B, type)
% DIST_TS_LOG: Measure the dissimilarity between two stochastic processes
%              based on covariance structures
%
% INPUT:
%   A: [vector] first time series
%   B: [vector] second time series
%   type: [string] 'cov': covariance matrix based distance
%                  'corr':  correlation matrix based distance
%
% OUTOUT:
%   d: distance between time series A and B

    if min(size(A,1), size(A,2)) ~= 1  
        error(' The first input series does not have dimension 1. \n');
    elseif max(size(A,1), size(A,2)) == 1
        error(' The first input series is a scalar, but should be a series. \n');
    elseif min(size(B,1), size(B,2)) ~= 1  
        error(' The second input series does not have dimension 1. \n');
    elseif max(size(B,1), size(B,2)) == 1
        error(' The first input series is a scalar, but should be a series. \n');
    end
    
    if length(A) >= length(B)
        p1 = A((end-length(B)+1):end); 
        p2 = B;
    else
        p1 = B((end-length(A)+1):end); 
        p2 = A;
    end
    
    if size(p1,1) ~= 1
        p1 = p1';
    end
    
    if size(p2,1) ~= 1
        p2 = p2';
    end
    
    n = length(p2);
    d = 0;
    
    for m = 1:floor(max(log(n),1));
        for L = 1:(n-m+1)
            nu_p1 = 0;
            nu_p2 = 0;
            for k = L:(n-m+1)
                nu_p1 = nu_p1 + 1/(n-m-L+2) * p1(k:(k+m-1))' * ...
                    p1(k:(k+m-1));
                nu_p2 = nu_p2 + 1/(n-m-L+2) * p2(k:(k+m-1))' * ...
                    p2(k:(k+m-1));
            end
            
            if strcmp(type, 'corr')
                if (trace(nu_p1) == 0) && (trace(nu_p2) == 0)
                    rho = 0;
                elseif ((trace(nu_p1) == 0) && (trace(nu_p2) ~= 0)) ...
                        || ((trace(nu_p1) ~= 0) && (trace(nu_p2) == 0))
                    rho = 1;
                else
                    % avoid numerical instability
                    rho = sqrt(max(1 - trace(nu_p1 * nu_p2) / ...
                        (sqrt(trace(nu_p1 * nu_p1)) * sqrt(trace(nu_p2 * nu_p2))),0));
                end
            elseif strcmp(type, 'cov')    
                rho = norm((sign_log(nu_p1) - sign_log(nu_p2)), 'fro');
            else
                rho = 0;
            end
            d = d + 1/(m*(m+1)) * 1/(L*(L+1)) * rho;
        end
    end
end


function log_sign_m=sign_log(m)
    log_sign_m = zeros(size(m));
    for i = 1:size(m, 1)
        for j = 1:size(m, 3)
            if m(i, j) > 0
                log_sign_m(i, j) = log(m(i, j));
            elseif m(i, j) < 0
                log_sign_m(i, j) = log(-m(i, j));
            end
        end
    end
end