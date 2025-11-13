function [varyingColsMatrix, varyingColsIndices] = findVaryingAttributesMatrix(X, cset, f)
    % X: The choice design matrix (nruns x nAttributes)
    % cset: Number of choice sets in the design
    % f: Number of fixed attributes in each choice set
    
    % Initialize parameters
    [nruns, nAttributes] = size(X);
    runsPerSet = nruns / cset;
    nVarying = nAttributes - f;
    
    % Initialize the vector to store varying columns indices
    varyingColsIndices = [];
    
    % Iterate through each choice set
    for csidx = 1:cset
        % Determine the start and end indices for the current choice set
        startIdx = (csidx - 1) * runsPerSet + 1;
        endIdx = csidx * runsPerSet;
        choiceSet = X(startIdx:endIdx, :);
        
        % Determine which attributes are varying (non-fixed)
        varyingCols = [];
        for colIdx = 1:nAttributes
            if length(unique(choiceSet(:, colIdx))) > 1
                varyingCols = [varyingCols, colIdx];
            end
        end
        
        % Store varying column indices
        varyingColsIndices = [varyingColsIndices; varyingCols'];
    end
    
    % Generate the varying columns matrix
    varyingColsMatrix = generateIndicatorMatrix(varyingColsIndices, nAttributes);
end


