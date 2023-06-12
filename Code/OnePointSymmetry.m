function [XX] = OnePointSymmetry(X,G,edges,d,i)
    
    [val1,ind1]=sort((edges(:,1)-i).^2);
    [val2,ind2]=sort((edges(:,2)-i).^2);
    ind1 = edges(ind1,2);
    ind2 = edges(ind2,1);

    indices = zeros(d,1);
    jj = 1;
    for j=1:length(ind1)
        if val1(j) == 0
            indices(jj) =ind1(j);
            jj = jj + 1;
        end

        if val2(j) == 0
            indices(jj) =ind2(j);
            jj = jj + 1;
        end
       
    end
    
    if degree(G,i) == d
        
        [Q,R] = qr((X(indices(2:size(indices,1)),:)-X(indices(1),:))');
        size(Q)
        coef_proj = X(i,:)*Q;
        x_star = X(i,:)';
        for k=1:size(Q,1)
            x_star = x_star - coef_proj(k)*Q(:,k);
        end
        x_star = -x_star
        for k=1:size(Q,1)
            x_star = x_star + coef_proj(k)*Q(:,k);
        end
       
    end

    XX = X;
    XX(i,:) = x_star';
end
