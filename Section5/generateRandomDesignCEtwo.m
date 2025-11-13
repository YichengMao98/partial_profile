function result = generateRandomDesignCEtwo(cset, n_alt, nlevels, varyingColsIndices)
    % Initialize parameters
    nAttributes = length(nlevels); % Number of attributes
    nVarying = length(varyingColsIndices) / cset; % Number of varying attributes per choice set

    % Initialize an empty cell array to hold data
    data = cell(1, nAttributes);
    
    % Loop through each cset
    for cs = 1:cset
        % Determine the start and end indices for the current choice set
        startIdx = (cs - 1) * nVarying + 1;
        endIdx = cs * nVarying;
        
        % Get the varying column indices for the current choice set
        currentVaryingCols = varyingColsIndices(startIdx:endIdx);

        % Determine the fixed columns as the complement of the varying columns
        fixedCols = setdiff(1:nAttributes, currentVaryingCols);
        
        % Loop through the levels
        for i = 1:nAttributes
            if ismember(i, fixedCols)
                % For fixed variables, assign a constant random value within their level range
                fixed_value = randi([1, nlevels(i)], 1, 1);
                data{i}((cs-1)*n_alt+1:cs*n_alt, 1) = repmat(fixed_value, n_alt, 1);
            else
                % For other variables, generate random data ensuring not all values are the same
                all_same = true;
                while all_same
                    random_values = randi([1, nlevels(i)], n_alt, 1);
                    if length(unique(random_values)) > 1
                        all_same = false;
                    end
                end
                data{i}((cs-1)*n_alt+1:cs*n_alt, 1) = random_values;
            end
        end
    end
    
    % Convert the cell array to a matrix
    result = cell2mat(data);
end