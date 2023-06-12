function [X, G, C, D, edges, A] = graphSphere3D(n,nn,m,random,sort_)
    if random
        if sort_
            theta = sort(unifrnd(0,2*pi,n,1));
            phi = sort(unifrnd(-pi/2,pi/2,nn+1,1));
        else
            theta = unifrnd(0,2*pi,n,1);
            phi = unifrnd(-pi/2,pi/2,nn+1,1);
        end
    else
        theta = 0:(2*pi/(n)):2*pi;
        theta = theta(1:n);

        phi = -pi/2:pi/(nn+2):pi/2;
        phi = phi(2:nn+2);
    end
    
    X = zeros(length(phi)*length(theta),3);
    for i=1:length(phi)
        phase = unifrnd(0,2*pi,1);
        for j=1:length(theta)
            X(j+(i-1)*length(theta),1) = cos(phi(i))*cos(theta(j)+phase);
            X(j+(i-1)*length(theta),2) = cos(phi(i))*sin(theta(j)+phase);
            X(j+(i-1)*length(theta),3) = sin(phi(i));
        end
    end
   
    if mod(m,2)>0
        m = m-1;
    end

    A = zeros(n,n);
    B = zeros(n,n);
    A(1,[2:1+m/2,1+n-m/2:n]) = 1;
    B(1,[1:1+m/2,1+n-m/2:n]) = 1;
    
    for i=2:n
        A(i,:) = circshift(A(i-1,:),1);
        B(i,:) = circshift(B(i-1,:),1);
    end

    A = kron(eye(nn+1),A);
    B = kron(eye(nn+1),B);
    A = triu(A + circshift(B,n,2)); % really bad behaviour with triu ... better without.
    A(n*(nn+1),1) = 1;
    A = A + A';
    A =  (A>0);
    G = graph(A);
    C = incidence(G);
    edges = table2array(G.Edges);
    D = eye(size(edges,1));
    for i=1:size(edges,1)
        D(i,i) = norm(X(edges(i,1),:)-X(edges(i,2),:),2);
    end
end