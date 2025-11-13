rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
num_ids = 30;

% #profiles
n_alt = 2;
% # fixed attribute
num_f = 2;
%  interactions
interactions = {};

% Define lambda and kappa values
lambda_values = [1, 1/2, 1/3];
kappa_values = [1, 1/2, 1/3];

% Preallocate results matrix (9 rows for 9 combinations)
results = zeros(9, 7);
index = 1; % To track row index in results

for lambda = lambda_values
    for kappa = kappa_values

        % Initialize design storage
        X_start = zeros(num_ids, cset*n_alt, length(nlevels));
        for id = 1:num_ids
            X = random_start_partial(cset, n_alt, nlevels, num_f);
            X_start(id, :, :) = X;
        end

        % Compute prior mean and variance
        priorMean = lambda * beta;
        priorVariance = (kappa^2) * sigma;
        [pts, wts] = gen_ptwt(priorMean, priorVariance);

        % Preallocate arrays for D_CE and time_CE
        D_best = -Inf;  % Initialize with negative infinity
        CE_best = [];   % Store best CE design
        time_CE_array = zeros(num_ids, 1);

        % Compute CE algorithm results
        for id = 1:num_ids
            start_time = tic;  

            X = squeeze(X_start(id, :, :));
            [~, varyingColsIndices] = findVaryingAttributesMatrix(X, cset, num_f);
            [varying_best, ~] = CE_stageone(varyingColsIndices, nlevels, cset, num_f);
            X_CE = generateRandomDesignCEtwo(cset, n_alt, nlevels, varying_best);
            [best_X_CE, best_D] = CE_stagetwo(X_CE, cset, n_alt, nlevels, num_f, interactions, pts, wts, varyingColsIndices);
            
            % Store best design in real-time
            if best_D > D_best
                D_best = best_D;
                CE_best = best_X_CE;  % Update best CE design
            end

            time_CE_array(id) = toc(start_time);
            fprintf('Completed processing for starting point %g in %g seconds\n', id, time_CE_array(id));
        end

        sum_time_CE = sum(time_CE_array); 

        % SA Algorithm Processing
        X = squeeze(X_start(1, :, :));
        X_main = transform_X(X, nlevels);
        D0 = calcObjFun_Bayesian(X_main, pts, wts, cset);
        delta = randomwalk_main(X, D0, cset, n_alt, nlevels, pts, wts, num_f, 100);
        T0 = abs(max(delta) / log(0.99));
        [best_X_SA, D_SA] = Bayesian_D_optimal_main_SA(X, cset, n_alt, nlevels, num_f, pts, wts, T0, sum_time_CE);

        % Compute Bayesian D-efficiency
        DB_eff = exp((D_best - D_SA) / length(priorMean));

        % Store results
        prior_int = length(priorMean) - 9; % Calculate interaction terms
        results(index, :) = [n_alt, num_f, prior_int, lambda, kappa, DB_eff, sum_time_CE];
        disp(results);
        
        % Save the best CE and SA designs with updated filenames
        ce_filename = sprintf('best_CE_alt%d_f%d_int%d_l%.2f_k%.2f.csv', n_alt, num_f, prior_int, lambda, kappa);
        sa_filename = sprintf('best_SA_alt%d_f%d_int%d_l%.2f_k%.2f.csv', n_alt, num_f, prior_int, lambda, kappa);
        writematrix(CE_best, ce_filename);
        writematrix(best_X_SA, sa_filename);
        
        fprintf('Saved best CE design to %s\n', ce_filename);
        fprintf('Saved best SA design to %s\n', sa_filename);
        
        index = index + 1; % Increment row index
    end
end

% Save results to CSV file
writematrix(results, 'alt2_f2_int0_results.csv');  

