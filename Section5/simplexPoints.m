function  simplexp =simplexPoints(p)
swt = ones(2 * p + 2,1) * (p * (7 - p)) / (2 * (p + 1) ^ 2 * (p + 2));
mwt = ones((p + 1) * p,1)* (2 * (p - 1) ^ 2) / (p * (p + 1) ^ 2 * (p + 2));
wt = [swt;mwt];
v = zeros(p + 1, p);
for i = 1:(p + 1)
        for j = 1:p
            if j < i
                v(i, j) = -sqrt((p + 1) / (p * (p - j + 2) * (p - j + 1)));
            elseif j == i
                v(i, j) = sqrt(((p + 1) * (p - i + 1)) / (p * (p - i + 2)));
            else
                v(i, j) = 0;
            end
        end
end
m = zeros(((p + 1) * p) / 2, p);
cCount = 1;
    for i = 1:p
        for j = (i + 1):(p + 1)
            tmp = (v(i, :) + v(j, :)) / 2;
            d = sqrt(tmp * tmp');
            m(cCount, :) = tmp / d;
            cCount = cCount + 1;
        end
    end
pts = [v;-v; m;-m];
simplexp=[pts,wt];
end