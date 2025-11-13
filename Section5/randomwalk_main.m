function delta = randomwalk_main(X,D0,cset,n_alt,nlevels,pts, wts,f,N)
current_X = X;
current_D = D0;
delta=[];
t = 1;
 while t <= N
       [row_idx,X_new] = random_modify_main(current_X, n_alt,nlevels,f);
       X_new_code = transform_X(X_new,nlevels);
       new_D = calcObjFun_Bayesian(X_new_code, pts, wts, cset);
       delta(:,t) = abs(new_D-current_D);
       t = t + 1;
       current_X = X_new;
       current_D = new_D;

 end
 delta = delta(delta <= 100);

end
