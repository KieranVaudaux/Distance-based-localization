%% Pathological case of local minima
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 20;
mm = 4 ;
r_minus_d = 0;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/NoRand_BadLocalMinima_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/NoRand_BadLocalMinima_d' ...
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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);

%% Fig 2
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 20;
mm = 4 ;
r_minus_d = 0;
fig = 2;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 3

legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 20;
mm = 4 ;
r_minus_d = 0;
fig = 3;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 4
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 50;
mm = 8 ;
r_minus_d = 0;
fig = 4;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 5
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 300;
mm = 8 ;
r_minus_d = 0;
fig = 5;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/NoRand_BadLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 5
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 50;
mm = 8 ;
r_minus_d = 0;
fig = 5;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/Rand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Rand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphCircle2D(n,mm,true,true);

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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 6
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 50;
mm = 14 ;
r_minus_d = 0;
fig = 6;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/Rand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Rand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'_rmd_',int2str(r_minus_d),'.png');
end

[X, G, C, D, edges, A] = graphCircle2D(n,mm,true,true);

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
 
title_ = strcat('Circular graph in 2D: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 7
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 20;
mm = 5 ;
r_minus_d = 0;
fig = 7;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/Rand_BadLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Rand_BadLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Graph on the disk with local structure: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
%% Fig 8
legend1 = 'Ground truth ';
legend2 = 'Found solution';
save = true;

d = 2;
n = 30;
mm = 5;
r_minus_d = 0;
fig = 8;

if r_minus_d ==0
    name = strcat('FigData/PathologicalCase/NoRelax/Rand_BadNoLocalMinima',int2str(fig),'_d' ...
        ,int2str(d),'_n_',int2str(n),'_m_',int2str(mm),'.png');
else
    name = strcat('FigData/PathologicalCase/Relax/Rand_BadNoLocalMinima',int2str(fig),'_d' ...
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
 
title_ = strcat('Graph on the disk with local structure: cost = ',num2str(xcost,3));
figure = plotGraph(G,X,x*x_rot, legend1, legend2, title_, save, name);
