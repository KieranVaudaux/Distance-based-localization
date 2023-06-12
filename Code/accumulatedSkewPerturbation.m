function [x_new,yy] = accumulatedSkewPerturbation(G,C,XX,i,num_pert)
% Method which try accumulate random perturbation for the vertex i and
% these neighbours

    N = neighbors(G,i);
    x_new = XX;
    for  j=1:num_pert
        [x_new,yy] = SphereIntersection(G,C,x_new,i,true);
        
        for l = 1:length(N)
            [x_new,yy] = SphereIntersection(G,C,x_new,l,true);
        end
    end

end