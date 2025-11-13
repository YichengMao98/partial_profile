function info = InfoMNL(X, b, cset)
    [r, c] = size(X);
    nchoices = r / cset;
    info = zeros(c, c);

    for i = 0:cset-1
        xrow = X(i * nchoices + 1 : (i + 1) * nchoices, :);
        u = xrow * b;
        exp_u = exp(u);
        p = exp_u ./ sum(exp_u);
        Pdiag = diag(p);
        info = info + xrow' * (Pdiag - p * p') * xrow;
    end
end
