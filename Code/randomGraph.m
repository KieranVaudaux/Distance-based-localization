function [X, G, C, D, edges, adj] = randomGraph(X,ratio)
    [n,d] = size(X);
    G = (rand(n)<ratio);
    G = (G+G'>0);
    G = graph(G,'omitselfloops');
    C = incidence(G);
    adj = adjacency(G);
    edges = table2array(G.Edges);
    D = eye(size(edges,1));

    for i=1:size(edges,1)
        D(i,i) = norm(X(edges(i,1),:)-X(edges(i,2),:),2);
    end

end