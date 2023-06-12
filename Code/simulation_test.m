%% Initialization of the problem

d = 2;
n = 50;
nn = 4;
mm = 15 ;

R=1;
r = 0.5;

X = normrnd(0,1,n,d);
[X, G, C, D, edges, A] = graphSphere(n,mm,d,true,false);%graphTorus3D(n,nn,mm,R,r,false,true);%graphCircle2D(n,mm,true,false);%graphSphere(n,mm,d,true);%graphCircle2D(n,mm,true,false);%graphSphere(n,mm,d,true);%graphTorus3D(n,nn,mm,R,r,false,true);%graphSphere3D(n,nn,mm,false,true);% randomGraph(X,ratio);graphSphere(n,mm,d,true);%
%[X, G, C, D, edges, A] = duplicateGraphHigherDimension(G,X,n,d,r_minus_d);

m = size(D,1)
%r = random('Normal',0,0.01,size(D,1));
D = D ;%+ diag(diag(r));
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;

% Define the problem cost function and its Euclidean gradient.
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
egrad = @(Y_) 2*Q2*Y_;   

S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));

%%
r_minus_d = 3;

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;
option.maxtime = 30;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

error = false;
[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,error,r);
[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj)

% x_new = x;
% for j=1:1000
%     i = randi(size(X,1));
%     [x_new,Y_perturbate] = SphereIntersection(G,C,x_new,i,false);
%     norm(x_new-x)
% end

%[x_new,Y_perturbate] = accumulatedSkewPerturbation(G,C,x,1,10);
%norm__ = norm(x-x_new)
%[x_proj_perturbate, diff_perturbation_perturbate, diff_perturbation_proj_perturbate] = deviationFromMeasureDistance(d,C,W,D,edges,Y_perturbate,error,r);


[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj);


figure = plotGraph(G,X,x_proj*x_rot,'','');
%figure = plotGraph(G,x_proj_perturbate,x_proj);
% diff_perturbation_perturbate, diff_perturbation_proj_perturbate

% cost_proj = cost(normr(C'*x_proj_perturbate))
% true_cost = cost(Y)

S = @(Y) Q2 - diag(diag(Q2*Y*Y'));
SS = S(Y);
min_eigs = min(eigs(SS,m))

%%
sigma = [];
iter = 1;
k=1;


for j=1:size(info_path_optimization,2)
    if isfield(info_path_optimization{k},"grad")
        sigma(iter,1) = info_path_optimization{k}.singular_val(1);
%         iter = iter + 1;
    end
    k = k + 1;
end

fig = semilogy(sigma,'DisplayName',int2str(1));
hold on;

for l=2:d+r_minus_d
    sigma = [];
    iter = 1;
    k=1;
    
    for j=1:size(info_path_optimization,2)
        if isfield(info_path_optimization{k},"grad")
            sigma(iter,1) = info_path_optimization{k}.singular_val(l);
            iter = iter + 1;
        end
        k = k + 1;
    end
    semilogy(sigma,'DisplayName',int2str(l));
    hold on;
end

semilogy(norm(Y_0)*2^(mean(degree(G))/n/log(d+1)).^(-(1:iter)),'--')
%semilogy(norm(Y_0)*2^(min(degree(G))/n/log(d+1)).^(-(1:iter)),'--')
%semilogy(norm(Y_0)*2^(0.5*log(min(degree(G)/d))*min(degree(G))/n).^(-(1:iter)),'--')
legend()
hold off;

% %%
% cost_ = [];
% sim_ = info_path_optimization{length(info_path_optimization)}.right_singular_space;
% iter = 1;
% k=1;
% 
% for j=1:size(info_path_optimization,2)
%     if isfield(info_path_optimization{k},"grad")
%         cost_(iter,1) = info_path_optimization{k}.cost_proj;
%         iter = iter + 1;
%     end
%     k = k + 1;
% end
% semilogy(cost_);
% %%
% r_minus_d =0;
% init = true;
% manifold = obliquefactory(d+r_minus_d, m,true);
% Z = 0.0001*ones(size(X_0,1),size(X_0,2)+1);
% Z(:,1:size(X_0,2)) = Y;
% Z = manifold.retr(Z, Z, 0);
% %X_0 = manifold.rand();
% option.tolgradnorm = 10^-9;
% 
% [x_up, Y_up, xcost_up, info_up, options_up] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Z,option);
% 
% error = false;
% [x_proj_up, diff_perturbation_up, diff_perturbation_proj_up] = deviationFromMeasureDistance(d,C,W,D,edges,Y_up,error,r);
% 
% [x_rot_up, xcost_rot_up, info_rot_up, options_rot_up] = rotationOptimization(d,X,x_proj_up);
% 
% %figure = plotGraph(G,X,x_proj_up*x_rot_up);
% 
% diff_perturbation_up, diff_perturbation_proj_up
% 
% S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));
% SS_up = S(Y_up);
% min(eigs(SS_up,m))
% 
% 
