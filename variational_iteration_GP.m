
function param = variational_iteration_GP(param_k, hyper, y, x)
%% Variational (blocked) coordinate ascent. 
%q() is factorized, we update each factor in turn. 
param = param_k; 
N = length(param.f.mu); 
 
 
 
%% (a) Recompute q_F(F)


P = hyper.kernel.K(x, hyper.alpha, N);
L = chol(P+hyper.fudge*eye(size(P)), 'upper')'; 
C = eye(N)/L; % L = C^{-1}


sqrtD = sqrt(param.stats.EZ/hyper.sigma2); 
R = triu(qr([bsxfun(@times, sqrtD, [L, y]) ; 
             eye(N), zeros(N, 1)])); 

param.f.cholprec = R(1:N, 1:N)*C; 

cholcovAndMu = param.f.cholprec\[eye(N), R(1:N, N+1)];
param.f.cholcov = cholcovAndMu(1:N, 1:N); 
param.f.mu = cholcovAndMu(1:N, N+1); 
param.f.cov = param.f.cholcov * param.f.cholcov'; 

assert(all(diag(param.f.cov)>0));
assert(all(~isnan(param.f.mu))); 
assert(all(isreal(param.f.mu))); 
assert(all(isreal(param.f.cholcov(:)))); 


param = updateParamStatsF(param); 



%% (b) Recompute q_Z(Z)
if hyper.type == 't'
    param.z.a = (hyper.nu + 1)/2;

    param.z.b = hyper.nu/2 ...
        + 1/(2*hyper.sigma2) * ((y-param.f.mu).^2 + diag(param.f.cov));

elseif hyper.type == 'b' 
    param.z.a = hyper.a + 1/2;
    param.z.b = hyper.b;
    param.z.k = ...
        - 1/(2*hyper.sigma2) * ((y-param.f.mu).^2 + diag(param.f.cov));
elseif hyper.type == 'g'
    % no update
    param.z.b = zeros(N, 1); 
    
else 
    error('unknown noise model!'); 
end

param = updateParamStatsZ(param, hyper.type); 




end