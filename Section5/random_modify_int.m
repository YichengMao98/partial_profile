function [csRows, modified_X] = random_modify_int(X, n_alt,nlevels,f,interactions)
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
        % Create set A with all interactions involving col_idx
           A = {};
           for i = 1:length(interactions)
               if ismember(col_idx, interactions{i})
                  A{end+1} = interactions{i};  % Add to set A
               end
           end
        % Create set B with all possible interactions between col_idx and other fixed_cols
           B = {};
           for j = 1:length(fixed_cols)
               if fixed_cols(j) ~= col_idx
                  interaction = sort([col_idx, fixed_cols(j)]);
                  B{end+1} = interaction;  % Add to set B
               end
           end
        % Check if all elements of A are in B
           interaction_fulfilled = all(cellfun(@(x) any(cellfun(@(y) isequal(x, y), B)), A));
        if interaction_fulfilled
           %change the level of this fixed attribute does not matters
           modified_X(row_idx, col_idx) = new_value;  
           % Randomly select another not fixed column and change it to fixed
           random_not_fixed_col = not_fixed_cols(randi(length(nlevels)-f));
           new_fixed_value = randi(nlevels(random_not_fixed_col));
           modified_X(csRows, random_not_fixed_col) = new_fixed_value;
           %disp(modified_X(csRows,:)-X(csRows,:));
        else
           p = f/length(nlevels);
           if rand() <= p
              % change all rows in the selected_rows with p
              modified_X(csRows, col_idx) = new_value;
           else
              modified_X(row_idx, col_idx) = new_value;  
              % Randomly select another not fixed column and change it to fixed
              random_not_fixed_col = not_fixed_cols(randi(length(nlevels)-f));
              new_fixed_value = randi(nlevels(random_not_fixed_col));
              modified_X(csRows, random_not_fixed_col) = new_fixed_value;
              %disp(modified_X(csRows,:)-X(csRows,:));
           end
        end
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

