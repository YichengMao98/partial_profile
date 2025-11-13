function ptwt = twoSpherePointsAndWeights(p)
pt0 = zeros(1, p);
wt0 = 8 / ((p + 2) * (p + 4));
wt1 = simplexPoints( p );
wt2 = wt1;
pt1 = wt1(:,1:p);
wt1(:, 1 : p) = [];
pt1 = pt1*sqrt( p + 4 - sqrt( 2 * p + 8 ) );
wt1 = wt1*p * (p + 2) / ((p + 4) * (2 - sqrt( 2 * p + 8 )) ^ 2);
pt2 = wt2(:,1:p);

wt2(:, 1 : p) = [];
pt2 = pt2*sqrt( p + 4 + sqrt( 2 * p + 8 ) );
wt2 = wt2*p * (p + 2) / ((p + 4) * (2 + sqrt( 2 * p + 8 )) ^ 2);
wt = [wt0 ;wt1 ; wt2];
pt = [pt0 ; pt1; pt2];
ptwt = [pt,wt];
end
