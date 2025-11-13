rng(2025); 
prior;
nlevels = [2,2,2,3,3,3];
cset = 24;
n_alt = 2;
num_f = 1;
X1 = design_main;
X2 = design_int;
X3 = design_robust;

%% main
priorMean_main = beta;
priorVariance_main = sigma;

[pts_main, wts_main] = gen_ptwt(priorMean_main, priorVariance_main);

X1_main = transform_X(X1,nlevels);
X2_main = transform_X(X2,nlevels);
X3_main = transform_X(X3,nlevels);

D1_main = calcObjFun_Bayesian(X1_main, pts_main,wts_main, cset);
D2_main = calcObjFun_Bayesian(X2_main, pts_main,wts_main, cset);
D3_main = calcObjFun_Bayesian(X3_main, pts_main,wts_main, cset);

disp(D1_main);%8.6718
disp(D2_main);%7.5197
disp(D3_main);%8.1744
exp((D2_main-D1_main)/9)%0.8798
exp((D3_main-D1_main)/9)%0.9462




