function [global_best_X, global_best_D] = Bayesian_D_optimal_int_SA(X,cset, n_alt, nlevels,f,interactions, pts, wts, initial_temp,max_time)
    start_time = tic;
    current_X = X;
    current_X_code_main = transform_X(current_X, nlevels);
    current_X_code = transform_X_inter(current_X_code_main,nlevels,interactions);
    current_infomats = geninfomats(current_X_code,pts,wts,cset);
    current_D = calcObjFun_Bayesian(current_X_code, pts, wts, cset);
    global_best_X = current_X;
    global_best_D = current_D;
    total_iter = 0;
    cycle = 1;
    cur_time = toc(start_time);
    while cur_time < max_time
          temp = initial_temp;
          t = 1;
          no_accept = 0;
          while no_accept <= 1000 && cur_time < max_time
            total_iter = total_iter + 1;

            [csRows,X_new] = random_modify_int(current_X, n_alt,nlevels,f,interactions);
            current_X_code_main = transform_X(current_X, nlevels);
            current_X_code = transform_X_inter(current_X_code_main, nlevels,interactions);
            X_new_code_main = transform_X(X_new,nlevels);
            X_new_code = transform_X_inter(X_new_code_main, nlevels,interactions);
            X_sel = current_X_code(csRows,:);
            X_sel_new = X_new_code(csRows,:);
            %new_D = calcObjFun_Bayesian(X_new_code, pts, wts, cset);
            infomats_new = zeros(size(X_sel,2),size(X_sel,2),length(wts));
            new_DB = 0;
            for idx = 1: length(wts)
                 b =pts(:,idx);
                 info_old = current_infomats(:,:,idx);
                 info_sel = InfoMNL(X_sel,b,1);
                 info_sel_new = InfoMNL(X_sel_new,b,1);
                 info_new = info_old+info_sel_new-info_sel;
                 infomats_new(:,:,idx) = info_new;
                 if det(info_new) >= 0 
                    new_DB = new_DB + log(det(info_new)) * wts(idx);
                 else
                    new_DB = -10000;
                 end
            end
            if new_DB >= current_D
               current_X = X_new;
               current_D = new_DB;
               current_infomats = infomats_new;
               no_accept = 0;
            elseif rand() < exp((new_DB-current_D) / temp)
               current_X = X_new;
               current_D = new_DB;
               current_infomats = infomats_new;
               no_accept = 0;
            else
               no_accept = no_accept + 1;
            end
            
            if current_D >  global_best_D
               global_best_X = current_X;
               global_best_D = current_D;
               fprintf('cycle: %d,total_iter: %d,time: %f, DB_best: %f\n',cycle,total_iter, cur_time,global_best_D);
            end
            t = t + 1;
            temp = initial_temp / t;
            cur_time = toc(start_time);
          end
    
    cycle = cycle + 1;
    end

end

