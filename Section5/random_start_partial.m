function result = random_start_partial(cset, n_alt, nlevels, f)
    % Initialize an empty cell array to hold data
    data = cell(1, length(nlevels));
    
    % Loop through each cset
    for cs = 1:cset
        % Select f random variables to fix
        fixed_vars = randperm(length(nlevels), f);
        
        % Loop through the levels
        for i = 1:length(nlevels)
            if ismember(i, fixed_vars)
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
