function xrow = eff_code(X, nlevels)
    nf = length(nlevels);
    df = sum(nlevels) - nf;
    xrow = zeros(1, df);
    startidx = 1; 
    
    for i = 1:nf
        nl = nlevels(i);
        xtmp = zeros(1, nl-1);
        
        if X(i) < nl
            xtmp(X(i)) = 1;
        else
            xtmp(1:nl-1) = -1;
        end
        
        xrow(startidx:startidx+nl-2) = xtmp;
        startidx = startidx + nl - 1;
    end
end

