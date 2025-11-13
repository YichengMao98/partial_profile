function X_code = transform_X(X, nlevels)
    [numRows, ~] = size(X);
    X_code = zeros(numRows, sum(nlevels) - length(nlevels)); % Pre-allocate matrix

    for row = 1:numRows
        X_code(row, :) = eff_code(X(row, :), nlevels);
    end
end