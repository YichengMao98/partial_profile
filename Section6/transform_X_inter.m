function interactionModel = transform_X_inter(X, nlevels, interactions)
    % 初始化交互模型矩阵
    interactionMatrix = [];

    % 计算每个属性的起始编码列索引
    colIndex = [1, cumsum(nlevels - 1) + 1];
    colIndex = colIndex(1:end-1); % 最后一个属性的起始索引不需要，因为累积总和已提供

    % 遍历每一对交互项
    for i = 1:length(interactions)
        cols = interactions{i};
        col1 = cols(1);
        col2 = cols(2);

        % 获取这两个属性的起始编码列索引
        startCol1 = colIndex(col1);
        startCol2 = colIndex(col2);
        
        % 获取这两个属性的编码列数（除去参考级别）
        nLevelsCol1 = nlevels(col1) - 1;  % 最后一个级别作为参考，不参与编码
        nLevelsCol2 = nlevels(col2) - 1;  % 同上

        % 生成这两个属性的所有有效级别的交互
        for level1 = 0:(nLevelsCol1 - 1)
            for level2 = 0:(nLevelsCol2 - 1)
                % 创建交互项列
                interactionColumn = X(:, startCol1 + level1) .* X(:, startCol2 + level2);
                % 将新的交互项列添加到交互模型矩阵
                interactionMatrix = [interactionMatrix, interactionColumn];
            end
        end
    end
    
    % 将结果作为逻辑矩阵返回
    interactionModel = [X,interactionMatrix];
end
