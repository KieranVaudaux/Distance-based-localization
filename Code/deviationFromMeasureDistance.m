function [X_star_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,error,r)
% Compute X_star from the solution and compare the distance between
% vertices obtains

X_star = pinv(C*W*C')*C*W*D*Y;
DD = diag(diag(D));
for i=1:size(DD,1)
    DD(i,i) = norm(X_star(edges(i,1),:)-X_star(edges(i,2),:),2);
end

[U,S,V] = svd(Y);
Y_proj = U*S(:,1:d);
Y_proj = normr(Y_proj);
X_star_proj = pinv(C*W*C')*C*W*D*Y_proj;
DD_proj = diag(diag(D));
for i=1:size(DD_proj,1)
    DD_proj(i,i) = norm(X_star_proj(edges(i,1),:)-X_star_proj(edges(i,2),:),2);
end
if error
    diff_perturbation = sum(abs(diag(D-DD) - diag(r)));
    diff_perturbation_proj = sum(abs(diag(D-DD_proj) - diag(r)));
else
    diff_perturbation = sum(abs(diag(D-DD)));
    diff_perturbation_proj = sum(abs(diag(D-DD_proj)));
end
end