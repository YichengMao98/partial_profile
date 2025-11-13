rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
n_alt = 2;
num_f = 1;
X = random_start_partial(cset, n_alt, nlevels, num_f);

priorMean_main = beta;
priorVariance_main = sigma;

[pts_main, wts_main] = gen_ptwt(priorMean_main, priorVariance_main);

X_main = transform_X(X,nlevels);

D0 = calcObjFun_Bayesian(X_main, pts_main,wts_main, cset);
delta = randomwalk_main(X,D0,cset,n_alt,nlevels, pts_main,wts_main,num_f,100);
T0 = abs(max(delta)/log(0.99));
[design_main, D_main,time_main] = Bayesian_D_optimal_main_SA(X,cset, n_alt, nlevels,num_f, pts_main,wts_main,T0);






