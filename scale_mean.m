function output = scale_mean(mat, mean_value)
% SCALE_MEAN: Scale each row of the input matrix to have specified mean
%   
% INPUT
%   mat: [matrix] input matrix;
%   mean_value: [scale] the user defined mean value.
%
% OUTPUT
%   output: [matrix] the scaled matrix whose mean of each row is defined by
%           'mean_value' in the input
output = mat;
for i = 1:size(mat, 1)
    output(i, :) = mat(i, :) - mean(mat(i, :)) + mean_value;
end

end

