function [pts,wts] = gen_ptwt(priorMean,priorVariance)
p = length(priorMean);
pw = twoSpherePointsAndWeights(p);
npts = size(pw,1);
wts = pw(:, p + 1)';
pts = pw(:,1:p)';
pv = chol( priorVariance );
pts = pv' * pts;
for bidx = 1:npts
    pts(:, bidx) = pts(:, bidx) + priorMean';
end