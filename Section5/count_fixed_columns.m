function num_fixed_cols = count_fixed_columns(matrix)
    % Function to count the number of fixed columns in a matrix
    % A fixed column is a column where all elements are the same
    
    % Initialize the counter
    num_fixed_cols = 0;
    
    % Loop through each column
    for col = 1:size(matrix, 2)
        % Check if all values in the current column are the same
        if length(unique(matrix(:, col))) == 1
            num_fixed_cols = num_fixed_cols + 1;
        end
    end
end