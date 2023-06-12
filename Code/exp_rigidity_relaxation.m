d = 2;
n = 30;
mm = 7 ;

%[X, G, C, D, edges, A] = graphSphere(n,mm,d,true,false);
X = normrnd(0,1,n,d);
[X, G, C, D, edges, A] = randomGraph(X,0.15);

m = size(D,1);
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_'); 

r_minus_d = 10;

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
option.tolgradnorm = 10^-6;
option.maxtime = 200;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[x_proj, diff_perturbation, diff_perturbation_proj] = deviationFromMeasureDistance(d,C,W,D,edges,Y,false,0);

[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj)

x_new = x;
for j=1:10
    i = randi(size(X,1));
    [x_new,Y_perturbate] = SphereIntersection(G,C,x_new,i,true);
    norm(x_new-x)
end

%[x_new,Y_perturbate] = accumulatedSkewPerturbation(G,C,x,2,500);
norm__ = norm(x-x_new)
[x_proj_perturbate, diff_perturbation_perturbate, diff_perturbation_proj_perturbate] = deviationFromMeasureDistance(d,C,W,D,edges,Y_perturbate,false,0);


xcost
cost_new = cost(normr(C'*x_new))
cost_proj
cost_pert = cost(normr(C'*x_proj_perturbate))
[x_rot, xcost_rot, info_rot, options_rot] = rotationOptimization(d,X,x_proj_perturbate);

figure = plotGraph(G,X,x_proj_perturbate*x_rot,'','','',false,'');
%title_ = strcat('d = ',int2str(d),', n = ',int2str(n),', m = ',int2str(mm),', costProj = ',num2str(cost_proj,3))