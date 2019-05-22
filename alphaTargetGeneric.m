
function [val, grad, hess] = alphaTargetGeneric(hyper, EfSq, x, alpha,n)

if nargout > 2
    [K, Kder, Khess] = hyper.kernel.K(x, alpha); 
elseif nargout > 1
    [K, Kder] = hyper.kernel.K(x, alpha); 
else
    K = hyper.kernel.K(x, alpha); 
end

np = length(alpha); 
Kfudge = K+hyper.fudge*eye(size(K)); 
Kinv = inv(Kfudge); 
KinvEfSq = Kinv*EfSq; %#ok<MINV> % K^{-1} EthetaSq

val = trace(KinvEfSq) + 2*sum(log(diag(chol(Kfudge)))); % F(alpha) = Tr[K(alpha)^-1 EthetaSq] + log det K(alpha)


assert(isreal(val))
if nargout>1
    traceIndex = bsxfun(@plus, 1:(n+1):(n^2), (0:np-1)'*n^2); 
    KinvKder = mtimes3(Kinv, Kder); % K^-1 K'
    
    KinvKderKinvEfSq = mtimes3(KinvKder, KinvEfSq); 
    
    
    grad = -sum(KinvKderKinvEfSq(traceIndex), 2) + sum(KinvKder(traceIndex), 2); 
    assert(isreal(grad))
% = -tr(P^-1 P_j P^-1 EthetaSq) + tr(P^-1 P_j)
    
%     grad = tr_j - tr_inv_j; % dF/d alpha_j = Tr( R_j EthetaSq ) - Tr(R(alpha)^-1 R_j)
    
    
   
    if nargout>2
        KinvKhess = mtimes3(Kinv, Khess); 
        
        hess1 = squeeze(sum(sum(bsxfun(@times, KinvKhess, ...
            permute(eye(size(K)) - KinvEfSq, [2 1 4 3])), 1), 2));
        hess2 = squeeze(sum(sum(bsxfun(@times, KinvKder, ...
            permute(2*KinvKderKinvEfSq - KinvEfSq, [2 1 4 3])), 1), 2));
        hess = hess1 + hess2; 
        hess = 0.5*(hess+hess'); 
        assert(isreal(hess))
    end
    
end

end