rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
n_alt = 2;
num_f = 1;
interactions = {[1,2],[1,4]};
priorMean_main = beta;
priorVariance_main = sigma;

X1 = design_main;
X2 = design_int;
X3 = design_robust;
X1_main = transform_X(X1, nlevels);
X2_main = transform_X(X2, nlevels);
X3_main = transform_X(X3, nlevels);
X1_int = transform_X_inter(X1_main, nlevels, interactions);
X2_int = transform_X_inter(X2_main, nlevels, interactions);
X3_int = transform_X_inter(X3_main, nlevels, interactions);

[pts_base, wts] = gen_ptwt(priorMean_main, priorVariance_main);
lambda_values = [0.1, 0.2, 0.3];

results = zeros(length(lambda_values), 4); 
for i = 1:length(lambda_values)
    lambda = lambda_values(i);

    new_rows = repmat(lambda, 3, size(pts_base, 2));
    pts = [pts_base; new_rows];

    D1 = calcObjFun_Bayesian(X1_int, pts, wts, cset);
    D2 = calcObjFun_Bayesian(X2_int, pts, wts, cset);
    D3 = calcObjFun_Bayesian(X3_int, pts, wts, cset);

    X = random_start_partial(cset, n_alt, nlevels, num_f);
    X_main = transform_X(X, nlevels);
    X_int = transform_X_inter(X_main, nlevels, interactions);
    D0 = calcObjFun_Bayesian(X_int, pts, wts, cset);
    delta = randomwalk_int(X, D0, cset, n_alt, nlevels, pts, wts, num_f, interactions, 100);
    T0 = abs(max(delta) / log(0.99));
    [design, D, time] = Bayesian_D_optimal_int_SA(X, cset, n_alt, nlevels, num_f, interactions, pts, wts, T0);

    eff1 = exp((D1 - D) / size(pts,1));
    eff2 = exp((D2 - D) / size(pts,1));
    eff3 = exp((D3 - D) / size(pts,1));
    results(i, :) = [lambda, eff1, eff2, eff3];

    fprintf('Lambda: %.2f | Eff1: %.4f | Eff2: %.4f | Eff3: %.4f\n', lambda, eff1, eff2, eff3);
end


