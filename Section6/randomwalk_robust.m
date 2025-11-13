function delta = randomwalk_robust(X,D0,cset,n_alt,nlevels,pts_main, wts_main, pts_full, wts_full,f,interactions,N)
current_X = X;
current_D = D0;
delta=[];
t = 1;
 while t <= N
       [row_idx,X_new] = random_modify_int(current_X, n_alt,nlevels,f,interactions);
       X_new_code_main = transform_X(X_new,nlevels);
       X_new_code_full = transform_X_inter(X_new_code_main,nlevels,interactions);
       %weight_main = size(pts_main,1)/(size(pts_main,1)+size(pts_full,1));
       %weight_full = size(pts_full,1)/(size(pts_main,1)+size(pts_full,1));
       %new_D = calcObjFun_Bayesian(X_new_code_main, pts_main, wts_main, cset)*weight_main+ calcObjFun_Bayesian(X_new_code_full, pts_full, wts_full, cset)*weight_full;
       new_D = calcObjFun_Bayesian(X_new_code_main, pts_main, wts_main, cset)/size(pts_main,1)+ calcObjFun_Bayesian(X_new_code_full, pts_full, wts_full, cset)/size(pts_full,1);
       delta(:,t) = abs(new_D-current_D);
       t = t + 1;
       current_X = X_new;
       current_D = new_D;

 end
end
