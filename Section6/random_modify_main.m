function [csRows, modified_X] = random_modify_main(X, n_alt,nlevels,f)
    [numRows, numCols] = size(X);
    row_idx = randi(numRows); 
    col_idx = randi(numCols);
    csRows = selectRows(row_idx, n_alt);
 
    modified_X = X; % Create a copy of X

    current_value = modified_X(row_idx, col_idx);
    new_value = generateRandomValue(current_value, nlevels(col_idx));
    
    % Check if the variable at col_idx is fixed or non-fixed
    all_cols = 1:numCols;
    fixed_cols = find(all(X(csRows(1), :) == X(csRows, :)));
    not_fixed_cols = setdiff(all_cols, fixed_cols);
    if ismember(col_idx, fixed_cols)
         modified_X(row_idx, col_idx) = new_value;  
         % Randomly select another not fixed column and change it to fixed
         random_not_fixed_col = not_fixed_cols(randi(length(nlevels)-f));
         new_fixed_value = randi(nlevels(random_not_fixed_col));
         modified_X(csRows, random_not_fixed_col) = new_fixed_value;
         %disp(modified_X(csRows,:)-X(csRows,:));
    else
        % If not fixed, change only the selected row
        modified_X(row_idx, col_idx) = new_value;
        if length(unique(modified_X(csRows, col_idx))) == 1
           % Randomly select another fixed column and change one of its values
           random_fixed_col = fixed_cols(randi(f));
           fixed_value = modified_X(csRows(1), random_fixed_col);
           new_value = generateRandomValue(fixed_value, nlevels(random_fixed_col));
           modified_X(csRows(randi(n_alt)), random_fixed_col) = new_value;
           %disp(modified_X(csRows,:)-X(csRows,:));
        end
    end
end
