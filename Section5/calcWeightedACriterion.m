function wtA = calcWeightedACriterion(XMatrix, nlevels)
    % Calculate x'x and then its inverse
    xtxi = inv(XMatrix' * XMatrix);

    % Initialize the weighted A-criterion accumulator
    wtA = 0;

    % Loop through each factor
    for idx = 1:length(nlevels)
        wtA = wtA + (((nlevels(idx) - 1)^2) / (2 * nlevels(idx))) * xtxi(idx, idx);
    end

    % Return the calculated weighted A-criterion
end
