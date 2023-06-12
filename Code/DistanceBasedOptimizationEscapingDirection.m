function [x, Y, cost, info, options] = DistanceBasedOptimizationEscapingDirection(D,C,W,d,Y_0,tolgrad, tolcost,max_iter)

m = size(D,1);

Q2 = -D*(W*C'*pinv(C*(W*C'))*C*W)*D;
Q1 = D*W*D;
manifold = obliquefactory(d, m,true);
problem.M = manifold;
 
manifold_up = obliquefactory(d+1, m,true);
problem_up.M = manifold_up;

% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(Y) trace(Q1 + Q2*Y*Y');
problem.egrad = @(Y) 2*Q2*Y;      
problem.ehess = @(Y,U) 2*Q2*U;

problem_up.cost  = @(Y) trace(Q1 + Q2*Y*Y');
problem_up.egrad = @(Y) 2*Q2*Y;      
problem_up.ehess = @(Y,U) 2*Q2*U;

criterion = true;
iter = 0;
Y = Y_0;

options.maxiter = 100;
options.tolgradnorm = tolgrad;
options_up.maxiter = 20;
options_up.tolgradnorm = tolgrad/10;

% Solve.
while criterion
    
    [Y, cost, info, ~] = trustregions(problem,Y,options);
    iter = iter + 1;

    if cost<tolcost 
        criterion = false;
    end
    current_gradnorm = problem.M.norm(Y,problem.M.egrad2rgrad(Y, problem.egrad(Y)));
    if current_gradnorm < tolgrad
        if cost>tolcost
            "inner loop"
            max_inner_iter = 20;
            cost_factor = 10;
            options_up.tolgradnorm = current_gradnorm/10;
            options_up.tolcost = cost/cost_factor;
            [Y, cost, inner_iter] = EscapingLocalMinima(Y,d,m,problem_up,options_up,cost,max_inner_iter,cost_factor);
            iter = iter + inner_iter;
        end
    end

    if iter>max_iter 
        criterion = false;
    end


end

x = pinv(C*W*C')*C*W*D*Y;