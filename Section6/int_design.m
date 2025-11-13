rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
n_alt = 2;
num_f = 1;
X = random_start_partial(cset, n_alt, nlevels, num_f);
%%int
interactions = {[1,2],[1,4]};
priorMean_int = [beta,zeros(1,3)];
priorVariance_int = blkdiag(sigma,eye(3));

[pts_int, wts_int] = gen_ptwt(priorMean_int, priorVariance_int);

X_main = transform_X(X,nlevels);
X_int = transform_X_inter(X_main,nlevels,interactions);
D0 = calcObjFun_Bayesian(X_int, pts_int,wts_int, cset);
delta = randomwalk_int(X,D0,cset,n_alt,nlevels,pts_int,wts_int,num_f,interactions,100);
T0 = abs(max(delta)/log(0.99));
[design_int, D_int,time_int] = Bayesian_D_optimal_int_SA(X,cset, n_alt, nlevels,num_f, interactions,pts_int,wts_int,T0);