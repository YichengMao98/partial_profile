rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
n_alt = 2;
num_f = 1;

X = random_start_partial(cset, n_alt, nlevels, num_f);
%%main
priorMean_main = beta;
priorVariance_main = sigma;

[pts_main, wts_main] = gen_ptwt(priorMean_main, priorVariance_main);

%%int
interactions = {[1,2],[1,4]};
priorMean_int = [beta,zeros(1,3)];
priorVariance_int = blkdiag(sigma,eye(3));

[pts_int, wts_int] = gen_ptwt(priorMean_int, priorVariance_int);

%%robust
X_main = transform_X(X,nlevels);
X_int = transform_X_inter(X_main,nlevels,interactions);

D0_robust = calcObjFun_Bayesian(X_main, pts_main, wts_main, cset)/size(pts_main,1)+ calcObjFun_Bayesian(X_int, pts_int, wts_int, cset)/size(pts_int,1);
delta = randomwalk_robust(X,D0_robust,cset,n_alt,nlevels,pts_main, wts_main, pts_int,wts_int,num_f,interactions,100);
T0 = abs(max(delta)/log(0.99));
[design_robust, D_robust,time_robust] = Bayesian_D_optimal_robust_SA(X,cset, n_alt, nlevels,num_f,interactions, pts_main, wts_main, pts_int,wts_int,T0);
