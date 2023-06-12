d = 3;
n = 30;
mm = 15 ;

[X, G, C, D, edges, A] = graphSphere(n,mm,d,false,true);
% X = normrnd(0,1,n,d);
% [X, G, C, D, edges, A] = randomGraph(X,0.2);

m = size(D,1);
W = eye(m);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;
cost  = @(Y_) trace(Q1 + Q2*Y_*Y_'); 

r_minus_d = 3;

init = true;
manifold = obliquefactory(d+r_minus_d, m,true);
Y_0 = manifold.rand();
[U,S,V] = svd(Y_0);
S(d+1:d+r_minus_d,d+1:d+r_minus_d) = S(d+1:d+r_minus_d,d+1:d+r_minus_d)*10^-1;
Y_0 = normr(U*S);
option.tolgradnorm = 10^-6;
option.maxtime = 200;

[x, Y, xcost, info, info_path_optimization] = DistanceBasedOptimization(D,C,W,d,r_minus_d,init,Y_0,option);

[U,S,V] = svd(Y);
Y_proj = normr(U*S(:,1:d));
cost_proj = cost(Y_proj)

title_ = strcat('d = ',int2str(d),', n = ',int2str(n),', m = ',int2str(mm),', costProj = ',num2str(cost_proj,3))
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

semilogy(norm(Y_0)*2^(mean(degree(G))/n/log2(d+1)).^(-(1:iter)),'--')
%semilogy(norm(Y_0)*2^(min(degree(G))/n/log(d+1)).^(-(1:iter)),'--')
%semilogy(norm(Y_0)*2^(0.5*log(min(degree(G)/d))*mean(degree(G))/n).^(-(1:iter)),'--')
legend()
title(title_)
hold off;
%saveas(fig,strcat('FigData/SingularValue/LocalStructBall/d_',int2str(d),'_n_',int2str(n),'_m_',int2str(mm ),'.png'))
