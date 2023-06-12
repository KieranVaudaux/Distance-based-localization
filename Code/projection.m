function [Y_proj] = projection(Y,d)
    [U,S,V] = svd(Y);
    Y_proj = U*S(:,1:d);
    Y_proj = normr(Y_proj);
end