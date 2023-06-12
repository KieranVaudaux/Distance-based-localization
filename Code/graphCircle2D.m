function [X, G, C, D, edges, A] = graphCircle2D(n,m,random,sort_)
    if random
        if sort_
            theta = sort(unifrnd(0,2*pi,n,1));
        else
            theta = unifrnd(0,2*pi,n,1);
        end
    else
        theta = 0:(2*pi/(n)):2*pi;
        theta = theta(1:n);
    end
    
    X = zeros(n,2);
    X(:,1) = sin(theta)';
    X(:,2) = cos(theta)';

    if mod(m,2)>0
        m = m-1;
    end

    A = zeros(n,n);
    A(1,[2:1+m/2,1+n-m/2:n]) = 1;
    
    for i=2:n
        A(i,:) = circshift(A(i-1,:),1);
    end

    
    G = graph(A);
    C = incidence(G);
    edges = table2array(G.Edges);
    D = eye(size(edges,1));
    for i=1:size(edges,1)
        D(i,i) = norm(X(edges(i,1),:)-X(edges(i,2),:),2);
    end
end