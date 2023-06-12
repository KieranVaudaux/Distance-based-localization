function [x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,X_star_proj)
% Runs an optimization algorithm rotate a given solution in order to try to
% match the ground truth solution

manifold = stiefelfactory(d, d, 1);
problem_.M = manifold;

n = size(X,1);
XX=X-repmat(mean(X),n,1);
X_star_proj = X_star_proj - mean(X_star_proj);

% Define the problem cost function and its Euclidean gradient.
problem_.cost  = @(QQ) trace((XX-X_star_proj*QQ)'*(XX-X_star_proj*QQ));
problem_.egrad = @(QQ) -(2*X_star_proj'*XX) + 2*X_star_proj'*X_star_proj*QQ; 

x_0 = manifold.rand();
% Solve.
[x_rot, xcost_rot, info_rot, options_rot] = trustregions(problem_,x_0);

if xcost_rot>10^-6
    x_0(:,1) = -x_0(:,1);
    [x_rot, xcost_rot, info_rot, options_rot] = trustregions(problem_,x_0);
end

end