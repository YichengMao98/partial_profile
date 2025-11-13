function indicatorMatrix = generateIndicatorMatrix(indices, nAttributes)
    % indices: A vector of column indices where indicators should be set to 1
    % nAttributes: Total number of attributes (number of columns in the matrix)
    
    % Initialize the indicator matrix
    numRows = length(indices);
    indicatorMatrix = zeros(numRows, nAttributes);
    
    % Set the appropriate columns to 1
    for i = 1:numRows
        indicatorMatrix(i, indices(i)) = 1;
    end
end
