%% Pathological case of local minima
%% Fig1
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 30;
mm = 4 ;
r_minus_d = 0;
fig = 1;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphCircle2D(n,mm,false,true);

m = size(D,1)
D = D ;
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;

% Define the problem cost function and its Euclidean gradient.
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
egrad = @(Y_) 2*Q2*Y_;   

S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x);

SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);

%%
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

r_minus_d = 1;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,false,r);
[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj);
S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));
SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, r = ',int2str(d+r_minus_d),' m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x_proj*x_rot, legend1, legend2, title_, save, name);


%% Fig2
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 100;
mm = 6 ;
r_minus_d = 0;
fig = 2;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphCircle2D(n,mm,false,true);

m = size(D,1)
D = D ;
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;

% Define the problem cost function and its Euclidean gradient.
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
egrad = @(Y_) 2*Q2*Y_;   

S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x);

SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);

%%
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

r_minus_d = 1;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,false,r);
[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj);
S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));
SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, r = ',int2str(d+r_minus_d),' m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x_proj*x_rot, legend1, legend2, title_, save, name);

%% Fig3
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 30;
mm = 6 ;
r_minus_d = 0;
fig = 3;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphSphere(n,mm,d,true,true);

m = size(D,1)
D = D ;
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;

% Define the problem cost function and its Euclidean gradient.
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
egrad = @(Y_) 2*Q2*Y_;   

S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x);

SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);

%%
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

r_minus_d = 2;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,false,r);
[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj);
S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));
SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, r = ',int2str(d+r_minus_d),' m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x_proj*x_rot, legend1, legend2, title_, save, name);


%% Fig4
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 50;
mm = 4 ;
r_minus_d = 0;
fig = 3;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphSphere(n,mm,d,true,false);

m = size(D,1)
D = D ;
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;

% Define the problem cost function and its Euclidean gradient.
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_');
egrad = @(Y_) 2*Q2*Y_;   

S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x);

SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);

%%
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

r_minus_d = 3;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Fig',int2str(fig),'/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,false,r);
[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj);

[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj);
S = @(Y_) Q2 - diag(diag(Q2*Y_*Y_'));
SS = S(Y);
min_eigs = min(eigs(SS,m));
 
title_ = strcat('d = 2, r = ',int2str(d+r_minus_d),' m = ',int2str(mm),', n = ',int2str(n),', cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x_proj*x_rot, legend1, legend2, title_, save, name);