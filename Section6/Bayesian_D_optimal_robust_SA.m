function [global_best_X, global_best_D,total_time] = Bayesian_D_optimal_robust_SA(X,cset, n_alt, nlevels,f,interactions, pts_main, wts_main, pts_full, wts_full,initial_temp)
    start_time = tic;
    current_X = X;
    %main model
    current_X_code_main = transform_X(current_X, nlevels);
    current_infomats_main = geninfomats(current_X_code_main,pts_main,wts_main,cset);
    current_D_main = calcObjFun_Bayesian(current_X_code_main, pts_main, wts_main, cset);
    %full model
    current_X_code_full = transform_X_inter(current_X_code_main, nlevels,interactions);
    current_infomats_full = geninfomats(current_X_code_full,pts_full,wts_full,cset);
    current_D_full = calcObjFun_Bayesian(current_X_code_full, pts_full, wts_full, cset);

    current_D = current_D_main/size(pts_main,1)+current_D_full/size(pts_full,1);
    global_best_X = current_X;
    global_best_D = current_D;

    total_iter = 0;
    cycle = 1;
    improved_in_cycle = true;
    while improved_in_cycle
          temp = initial_temp;
          t = 1;
          no_accept = 0;
          improved_in_cycle = false; % Reset improvement flag at the start of each cycle
          while no_accept <= 1000
            total_iter = total_iter + 1;

            [csRows,X_new] = random_modify_int(current_X, n_alt,nlevels,f,interactions);
            %main model
            current_X_code_main = transform_X(current_X, nlevels);
            X_new_code_main = transform_X(X_new,nlevels);
            X_sel_main = current_X_code_main(csRows,:);
            X_sel_new_main = X_new_code_main(csRows,:);
            infomats_new_main = zeros(size(X_sel_main,2),size(X_sel_main,2),length(wts_main));
            new_DB_main = 0;
            for idx = 1: length(wts_main)
                 b =pts_main(:,idx);
                 info_old = current_infomats_main(:,:,idx);
                 info_sel = InfoMNL(X_sel_main,b,1);
                 info_sel_new = InfoMNL(X_sel_new_main,b,1);
                 info_new = info_old+info_sel_new-info_sel;
                 infomats_new_main(:,:,idx) = info_new;
                 if det(info_new) >= 0 
                    new_DB_main = new_DB_main + log(det(info_new)) * wts_main(idx);
                 else
                    new_DB_main = -10000;
                 end
            end
            %full model
            current_X_code_full = transform_X_inter(current_X_code_main, nlevels,interactions);
            X_new_code_full = transform_X_inter(X_new_code_main, nlevels,interactions);
            X_sel_full = current_X_code_full(csRows,:);
            X_sel_new_full = X_new_code_full(csRows,:);
            infomats_new_full = zeros(size(X_sel_full,2),size(X_sel_full,2),length(wts_full));
            new_DB_full = 0;
            for idx = 1: length(wts_full)
                b =pts_full(:,idx);
                info_old = current_infomats_full(:,:,idx);
                info_sel = InfoMNL(X_sel_full,b,1);
                info_sel_new = InfoMNL(X_sel_new_full,b,1);
                info_new = info_old+info_sel_new-info_sel;
                infomats_new_full(:,:,idx) = info_new;
                if det(info_new) >= 0 
                   new_DB_full = new_DB_full + log(det(info_new)) * wts_full(idx);
                else
                   new_DB_full = -10000;
                end
            end
            new_DB = new_DB_main/size(pts_main,1) + new_DB_full/size(pts_full,1);
            if new_DB >= current_D
               current_X = X_new;
               current_D = new_DB;
               current_infomats_main = infomats_new_main;
               current_infomats_full = infomats_new_full;
               no_accept = 0;
            elseif rand() < exp((new_DB-current_D) / temp)
               current_X = X_new;
               current_D = new_DB;
               current_infomats_main = infomats_new_main;
               current_infomats_full = infomats_new_full;
               no_accept = 0;
            else
               no_accept = no_accept + 1;
            end
            
            if current_D >  global_best_D
               global_best_X = current_X;
               global_best_D = current_D;
               improved_in_cycle = true; % Improvement occurred in this cycle
               current_time = toc(start_time);
               fprintf('cycle: %d,total_iter: %d,time: %f, DB_best: %f\n',cycle,total_iter, current_time,global_best_D);
            end
        t = t + 1;
        temp = initial_temp / t;
          end
    cycle = cycle + 1;
    current_time = toc(start_time);
    fprintf('cycle: %d,total_iter: %d,time: %f, DB_best: %f\n',cycle,total_iter, current_time,global_best_D);
    end
    total_time = toc(start_time);
end

