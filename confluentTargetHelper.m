
function [val, grad, hess] = confluentTargetHelper(ab,N,EsumlogZ, Esumlog1mZ)
%     logExpectedPriorZ = @(a,b,N) (a-1)*param.stats.EsumlogZ + (b-1)*param.stats.Esumlog1mZ ...
%         - N*betaln(a,b); %
val = -(ab(1)-1)*EsumlogZ + -(ab(2)-1)*Esumlog1mZ + N*(sum(gammaln(ab))-gammaln(sum(ab))); 
if nargout > 1
    
    grad = [-EsumlogZ+N*(psi(0,ab(1))-psi(0,sum(ab)));
            -Esumlog1mZ+N*(psi(0,ab(2))-psi(0,sum(ab)))];
    hess = N*([psi(1,ab(1)), 0; 
            0, psi(1,ab(2))] - psi(1,sum(ab))); 
    
end
end