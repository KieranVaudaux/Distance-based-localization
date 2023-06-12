function [x_new,yy] = SphereIntersection(G,C,XX,i,skew)
% Random perturbation on the intersection of spheres
    N = neighbors(G,i);
    
    D_  = diag((XX(N,:) - XX(i,:))*(XX(N,:)-XX(i,:))');
    
    A = XX(N,:)';
    a_m = A(:,end);
    A_tilde = A(:,1:end-1) - a_m;
    
    col_norm_A = diag(A_tilde'*A_tilde);
    
    [Q,R_] = qr(A_tilde);
   
    k = rank(A_tilde,10^-8);
    R_tilde = R_(1:k,:);
    c = -0.5*(D_(1:end-1) - D_(end) - col_norm_A);
    
    y = R_tilde'\c;

    n = size(A,1);
    n_k = n - k;

    x_new = XX;

    x_yz = zeros(n,1);
    x_yz(1:k) = y;
    t = D_(end)-y'*y;
    if t < 0
        t = 0;
    end

    if skew
        x_yz(k+1:n) = sqrt(t)*normc(abs(normrnd(0,1,n_k,1)));
    else
        x_yz(k+1:n) = sqrt(t)*normc(normrnd(0,1,n_k,1));
    end

    x_new_ = Q*x_yz+ a_m;

    x_new(i,:) = x_new_';
    yy = normr(C'*x_new);

    %norm(D_ - diag((x_new(N,:) - x_new(i,:))*(x_new(N,:)-x_new(i,:))'))


end