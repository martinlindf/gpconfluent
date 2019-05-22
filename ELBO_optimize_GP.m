
function [hyper, ELBOmax] = ELBO_optimize_GP(param, hyper_k, y, x)
hyper = hyper_k; 
persistent opt opt2; 

% compute some useful quantities. See ELBO()
N = length(param.f.mu);
n = length(param.f.mu); 



if isempty(opt2)
%     opt2 = optimoptions('fmincon',...
%         ...'Algorithm','trust-region-reflective', ...
%         'Algorithm','interior-point', ...
%         'SpecifyObjectiveGradient', false,...
%         ...'HessianFcn','objective', ...
%         'OptimalityTolerance', 1e-3, 'FunctionTolerance', 1e-3, 'Display', 'off'...
%         );
    opt2 = optimoptions('fmincon',...
        'Algorithm','trust-region-reflective', ...
        ...'Algorithm','interior-point', ...
        'SpecifyObjectiveGradient', true,...
        'HessianFcn','objective', ...
        'OptimalityTolerance', 1e-5, 'FunctionTolerance', 1e-5, 'Display', 'off'...
        );
    
end

% update sigma^2
hyper.sigma2 = param.stats.EZ' * ((y-param.f.mu).^2 + diag(param.f.cov))/N;

assert(hyper.sigma2 > 0); 
% update eta_z
if hyper.type == 't'
    % update nu; explicit approximate expression (by inspection)
    % accurate within 0.01 of the correct optimum. 
    ebar = (param.stats.EsumZ-param.stats.EsumlogZ)/N;
    assert(ebar>1); 
    hyper.nu = 1/(ebar-1) + .44*(1-2*atan(ebar/2.3)/pi); 
    assert(hyper.nu > 0); 
elseif hyper.type == 'b'
    targetC = @(ab)confluentTargetHelper(ab,N,param.stats.EsumlogZ, param.stats.Esumlog1mZ);
    
%     a = fminbnd(@(a)targetC([a,hyper.b]), 0.51, 20); 
%     if isscalar(a)
%         hyper.a = a;
%     end
%     b = fminbnd(@(b)targetC([hyper.a,b]), 0, 1); 
%     if isscalar(b)
%         hyper.b = b;
%     end
    
    ab = fmincon(targetC, [hyper.a; hyper.b],...
        [],[],[],[],[0.01; 0], [20;1],[],opt2); 
    if numel(ab)==2
        hyper.a = ab(1);
        hyper.b = ab(2);
%         hyper.a = 1; 
%         hyper.b = 0.2; 
    end
    assert(hyper.a>0);
    assert(hyper.b>0); 
elseif hyper.type == 'g'
    % eta_z is empty
else 
    error('unknown noise model!'); 
end

% optimize alpha. 

% if the target function is explicitly specified, pick it. Else use helper
assert(~any(isnan(param.stats.EfSq(:)))); 
target = @(alpha) alphaTargetGeneric(hyper, param.stats.EfSq, x, alpha, n);

if isempty(opt)
%     if isempty(hyper.kernel.alphaCon) 
%     opt = optimoptions('fmincon',...
%         'Algorithm','trust-region-reflective', ...
%         ...'Algorithm','interior-point', ...
%         'SpecifyObjectiveGradient', true, 'HessianFcn','objective', ...
%         'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'Display', 'off'...
%         );
%     else
    opt = optimoptions('fmincon',...
        ...'Algorithm','trust-region-reflective', ...
        'Algorithm','interior-point', ...
        ...'SpecifyObjectiveGradient', true, ...
        ...'HessianFcn','objective', ...
        'OptimalityTolerance', 1e-3, 'FunctionTolerance', 1e-3, 'Display', 'off'...
        );
%     end
end

% newalpha = [-2 5]; 
% try
% [val, grad, hess] = target(hyper.alpha); 
    newalpha = fmincon(target, hyper.alpha, [],[],[],[], ...
            hyper.kernel.alphaLB, hyper.kernel.alphaUB, ...
            hyper.kernel.alphaCon, opt);
% catch err
%     sin(1);
% end


hyper.alpha = newalpha; 
ELBOmax = ELBO_GP(param, hyper, y, x);

 

end