
function [K, Kder, Khess] = kernelSE(x, alpha,n) 
% x: n x nx
x = permute(x, [1 3 2]); 
xperm = permute(x, [2 1 3]); 
 % Delta_i,j = sum((x_{i,k} - x_{j,k}).^2, k)
Delta = sum((x-xperm).^2, 3); 

scaledDelta = -Delta*exp(alpha(2)); 
K = exp(alpha(1)) .* exp(scaledDelta); 
if nargout > 1
    dKdalpha1 = K;
    dKdalpha2 = K .* scaledDelta;
    Kder = cat(3, dKdalpha1, dKdalpha2); 
    if nargout > 2
        dKdalpha1dalpha1 = K;
        dKdalpha1dalpha2 = dKdalpha2; 
        dKdalpha2dalpha1 = dKdalpha2; 
        dKdalpha2dalpha2 = dKdalpha2 .* (1 + scaledDelta); 
        Khess = cat(4, cat(3, dKdalpha1dalpha1, dKdalpha1dalpha2), ...
            cat(3, dKdalpha2dalpha1, dKdalpha2dalpha2)); 
        
    end
end

end
