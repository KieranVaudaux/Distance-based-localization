function [X_new, G_new, C, D, edges, A_new] = duplicateGraphHigherDimension(G,X,n,d,r_minus_d)
    % Take  graph in a embedded in a given dimension and create a new graph
    % in a higher dimension from the given one by duplicate him, then
    % rotate the duplicate and connect both.
    manifold = stiefelfactory(d+r_minus_d, d+r_minus_d, 1);%rotationsfactory(d, 1);
    Q = manifold.rand();


    X_new = zeros(2*size(X,1),d+r_minus_d);
    X_new(1:size(X,1),1:d) = X;
    X_new((size(X,1)+1):2*size(X,1),1:d) = X;

    X_new((size(X,1)+1):2*size(X,1),:) = X_new((size(X,1)+1):2*size(X,1),:)*Q+ 5;

    A_new = zeros(2*n);
    A = full(adjacency(G));
    A_new(1:n,1:n) = A;
    A_new((n+1):2*n,(n+1):2*n)= A;
    
    D = zeros(n);
    for i=1:n
        for  j=1:n
            D(i,j) = norm(X_new(i,:)-X_new(i+n,:));
        end
    end
    D_ = D(:);
    [val,ind_] = sort(D_);
    
    community_connection = d + r_minus_d + 20;
    for k=1:community_connection
        if mod(ind_(k),n) == 0
            
            i = ind_(k)/n ;
            j = n;
        else
            
            i = (ind_(k) - mod(ind_(k),n))/n + 1;
            j = mod(ind_(k),n);
        end
        A_new(i,j+n) = 1;
        A_new(j+n,i) = 1;

        
    end

    G_new = graph(A_new);
    C = incidence(G_new);
    edges = table2array(G_new.Edges);
    D = eye(size(edges,1));
    for i=1:size(edges,1)
        D(i,i) = norm(X_new(edges(i,1),:)-X_new(edges(i,2),:),2);
    end

end