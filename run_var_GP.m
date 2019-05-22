function [elbo, param_k, hyper_k, iterMax, elbo_hist, param_hist, hyper_hist] = ...
    run_var_GP(y, x, param_0, hyper_0, iterMin, iterMax, iterVarMax, tol)

param_k = param_0;
hyper_k = hyper_0; 
N = length(param_k.f.mu); 

iter = 1;
hyper_hist = cell(iterMax, 1); 
param_hist = cell(iterMax, iterVarMax); 
elbo_hist = NaN(iterMax , iterVarMax + 1);
done = false;  
while (~done && iter<=iterMax)
    %% Run a full iteration. 
    
    % (1) Do a few variational iterations to fit the variational posterior
    % approximation q(f, z) = param to p(f, z|y,alpha)
    
    % (2) Do an optimization over the ELBO to estimate eta = hyper
    iterVar = 1; 
    while (iterVar<=iterVarMax)
        %% (1) Variational iterations
        param_k = variational_iteration_GP(param_k, hyper_k, y, x); 
        
        param_hist{iter, iterVar} = param_k; 
        elbo_hist(iter, iterVar) = ELBO_GP(param_k, hyper_k, y, x); 
        iterVar = iterVar + 1; 
    end
    
    %% (2) ELBO optimization
    
    [hyper_k, elbo] = ELBO_optimize_GP(param_k, hyper_k, y, x); 
    hyper_hist{iter} = hyper_k; 
    
    elbo_hist(iter, end) = elbo; 
    
    
    done = (iter>=iterMin) ...
            && ((elbo_hist(iter, end) - elbo_hist(iter-1, end))/N < tol); 
    
    iter = iter+1; 
end

iterMax = iter-1;
end