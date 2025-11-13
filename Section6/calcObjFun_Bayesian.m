function objfun = calcObjFun_Bayesian(xmat, pts, wts, cset)
    objfun = 0;
 
    for idx = 1:length(wts)
        b = pts(:, idx);
        info = InfoMNL(xmat, b, cset);
        det_info = det(info);
        
       if det_info >= 0
            objfun = objfun + log(det_info) * wts(idx);
       else
            objfun = -10000;
       end
       
    end
end
