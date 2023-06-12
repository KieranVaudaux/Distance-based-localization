function [X, G, C, D, edges, A] = graphTorus3D(n,m,R,r,local)
    theta = unifrnd(0,2*pi,n,1);
    phi = unifrnd(0,2*pi,n,1);
    
    X = zeros(n,3);
    for i=1:length(phi)
        X(i,1) = cos(phi(i))*(R+r*cos(theta(i)));
        X(i,2) = sin(phi(i))*(R+r*cos(theta(i)));
        X(i,3) = r*sin(theta(i));
    end

    D = zeros(n,n);
    for i=1:n
        for j=i:n
            D(i,j) = norm(X(i,:)-X(j,:),2);
        end
    end
    D = D + D';
    A = zeros(n,n);
    
    for i=1:n
        [sortval,ind] = sort(D(i,:));
        if not(local)
             A(i,ind(n-1)) = 1;
             A(ind(n-1),i) = 1;
        end
        if sum(A(i,:))<m
            [sortval,ind] = sort(D(i,:));
            for r=2:m
                A(i,ind(r)) = 1;
                A(ind(r),i) = 1;
            end
        end
    end



    G = graph(A);
    C = incidence(G);
    edges = table2array(G.Edges);
    D = eye(size(edges,1));
    for i=1:size(edges,1)
        D(i,i) = norm(X(edges(i,1),:)-X(edges(i,2),:),2);
    end
end