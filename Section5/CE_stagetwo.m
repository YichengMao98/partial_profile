function [global_best_X, global_best_D] = CE_stagetwo(X, cset, n_alt, nlevels, f, interactions, pts, wts, varyingColsIndices)
    nruns = cset * n_alt;
    nAttributes = length(nlevels);
    nVarying = nAttributes - f;
    current_X = X;
    current_X_code_main = transform_X(current_X, nlevels);
    current_X_code = transform_X_inter(current_X_code_main, nlevels, interactions);
    current_infomats = geninfomats(current_X_code, pts, wts, cset);
    current_D = calcObjFun_Bayesian(current_X_code, pts, wts, cset);
    global_best_X = current_X;
    global_best_D = current_D;

    madeswitch = 1;
    iter = 0;
    while madeswitch > 0 
        madeswitch = 0;
        for row = 1:nruns
            csidx = ceil(row / n_alt);
            cs_rows = (csidx - 1) * n_alt + 1 : csidx * n_alt;
            startIdx = (csidx - 1) * nVarying + 1;
            endIdx = csidx * nVarying;
            VaryingCols = varyingColsIndices(startIdx:endIdx);

            for col = 1:nAttributes
                new_X = current_X;
                levels = 1:nlevels(col);
                combinations = nchoosek(levels, n_alt);
                for comb_idx = 1:size(combinations, 1)
                    combination = combinations(comb_idx, :);
                    permutations = perms(combination);
                    for perm_idx = 1:size(permutations, 1)
                        iter = iter + 1;
                        perm = permutations(perm_idx, :);

                        if ismember(col, VaryingCols)
                           for k = 1:n_alt
                               new_X(cs_rows(k), col) = perm(k);
                           end
                
                           if length(unique(new_X(cs_rows, col))) == 1
                               continue;
                           end
                        else
                           new_X(cs_rows, col) = perm(1);
                        end
                        num_fixed_cols = count_fixed_columns(new_X(cs_rows, :));
                        if num_fixed_cols ~= f
                           continue;
                        end
                        current_X_code_main = transform_X(current_X, nlevels);
                        current_X_code = transform_X_inter(current_X_code_main, nlevels, interactions);
                        X_new_code_main = transform_X(new_X, nlevels);
                        X_new_code = transform_X_inter(X_new_code_main, nlevels, interactions);
                        X_sel = current_X_code((csidx - 1) * n_alt + 1 : csidx * n_alt, :);
                        X_sel_new = X_new_code((csidx - 1) * n_alt + 1 : csidx * n_alt, :);
                        new_DB = 0;
                        infomats_new = zeros(size(X_sel, 2), size(X_sel, 2), length(wts));
                        for idx = 1:length(wts)
                            b = pts(:, idx);
                            info_old = current_infomats(:, :, idx);
                            info_sel = InfoMNL(X_sel, b, 1);
                            info_sel_new = InfoMNL(X_sel_new, b, 1);
                            info_new = info_old + info_sel_new - info_sel;
                            infomats_new(:, :, idx) = info_new;
                            if det(info_new) >= 0
                               new_DB = new_DB + log(det(info_new)) * wts(idx);
                            else
                               new_DB = -10000;
                            end
                        end

                        if new_DB > current_D
                           current_X = new_X;
                           current_D = new_DB;
                           current_infomats = infomats_new;
                           madeswitch = 1;
                        end
                        if current_D > global_best_D
                           global_best_X = current_X;
                           global_best_D = current_D;
                        end
                   end
               end
           end
       end
    end
end



