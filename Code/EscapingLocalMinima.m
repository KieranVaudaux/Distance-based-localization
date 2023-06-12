function [Y_proj, inner_cost, inner_iter] = EscapingLocalMinima(Y,d,m,inner_problem,inner_options,current_cost,max_inner_iter,cost_factor)
% Algorithm whic increase the dimension of the problem when we converge to
% local minima, continue the optimization until be far from the previous
% old minima and then project back in the original space.
Y_ = zeros(m,d+1);
Y_(:,1:d) = Y;
Y_(:,d+1) = normrnd(0,inner_options.tolgradnorm*100,m,1);
Y_ = normr(Y_);
[Y_, cost_, info_, ~] = trustregions(inner_problem,Y_,inner_options);

Y_proj = projection(Y_,d);
inner_cost = inner_problem.cost(Y_proj);

inner_iter = 1;
inner_criterion = false;

while (inner_criterion & cost_>current_cost/cost_factor & inner_iter<max_inner_iter)
    [Y_, cost_, info_, ~] = trustregions(inner_problem,Y_,inner_options);
    
    inner_iter = inner_iter + 1
    cost_>current_cost/cost_factor
    if cost_>current_cost/cost_factor
        Y_proj = projection(Y_,d);
        inner_cost = inner_problem.cost(Y_proj);
        if inner_cost<current_cost
            inner_criterion = false;
        end
    end 
end


end