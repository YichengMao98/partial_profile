function [varying_best,wtA_best]= CE_stageone(varyingColsIndices,nlevels , cset, f)
    nAttributes =length(nlevels);
    nVarying = nAttributes - f;
    
    %%%Z
    repeatingVector = zeros(cset * nVarying, 1);

    % Fill the vector
    for i = 1:cset
        startIdx = (i - 1) * nVarying + 1;
        endIdx = i * nVarying;
        repeatingVector(startIdx:endIdx) = i;
    end
    Z = transform_X(repeatingVector,cset);
    ZtZ = Z' * Z;
    ZtZ_inv = inv(ZtZ);
    mid_matrix = eye(cset* nVarying)-Z * ZtZ_inv * Z';

    Q0 = generateIndicatorMatrix(varyingColsIndices, nAttributes);
    N_matrix0 = Q0' * mid_matrix * Q0;
    wtA0 = calcWeightedACriterion(N_matrix0, nlevels);
    wtA_best = wtA0;
    varying_best = varyingColsIndices;

    madeswitch = 1;
    iter = 0;
    while madeswitch > 0 
        madeswitch = 0;
        % Generate all combinations of df choose (df-f)
        combinations = nchoosek(1:nAttributes, nVarying);
        for csidx = 1:cset
        % Iterate through each combination
            for combIdx = 1:size(combinations, 1)
            % Get the current combination of varying columns
                iter = iter+1;
                combVaryingCols = combinations(combIdx, :);

            % Update the varyingColsIndices for the current choice set
                startIdx = (csidx - 1) * nVarying + 1;
                endIdx = csidx * nVarying;
                updatedVaryingColsIndices = varyingColsIndices;
                updatedVaryingColsIndices(startIdx:endIdx) = combVaryingCols;

            % Generate the indicator matrix for the updated varying columns indices
                Q = generateIndicatorMatrix(updatedVaryingColsIndices, nAttributes);
            
            % Compute the N matrix
               N_matrix = Q' * mid_matrix * Q;
               wtA = calcWeightedACriterion(N_matrix, nlevels);
               if wtA < wtA_best
                  wtA_best = wtA;
                  varying_best = updatedVaryingColsIndices;
                  madeswitch = 1;
               end
            end
        
        end
       
    end
end

