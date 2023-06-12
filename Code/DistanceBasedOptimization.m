function [x, Y, xcost, info,info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,X_0,opt)
% Function which runs a optimization for a given graph and set of measured
% distances

m = size(D,1);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;
manifold = obliquefactory(d+r_minus_d, m,true);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(Y) trace(Q1 + Q2*Y*Y');
problem.egrad = @(Y) 2*Q2*Y;      
problem.ehess = @(Y,U) 2*Q2*U;

% Solve.
if init
    [Y, xcost, info, options, info_path_optimization] = trustregions_trackSubSpaceDim(problem, X_0, opt,d,d+r_minus_d,problem.cost);
else
    [Y, xcost, info, options] = trustregions_trackSubSpaceDim(problem, problem.M.rand(), opt,d,d+r_minus_d,problem.cost);
end
x = pinv(C*W*C')*C*W*D*Y;
end