
function [K, Kder, Khess] = kernelSEnD(x, alpha,n) 
% x: n x nx
nx = size(x, 2); 
x = permute(x, [1 3 2]); 
xperm = permute(x, [2 1 3]); 
 % Delta_i,j = sum((x_{i,k} - x_{j,k}).^2, k)
Delta = ((x-xperm).^2); 
loginvscale = permute(alpha(2:nx+1), [3 2 1]); % one per dimension
 
scaledDelta = -Delta.*exp(loginvscale); 
K = exp(alpha(1)) .* exp(sum(scaledDelta, 3)); 
if nargout > 1
    dKdalpha1 = K;
    dKdalpha2 = K .* scaledDelta;
    Kder = cat(3, dKdalpha1, dKdalpha2); 
    if nargout > 2
        dKdalpha1dalpha1 = K;
        dKdalpha1dalpha2 = dKdalpha2; 
        dKdalpha2dalpha1 = permute(dKdalpha2, [1 2 4 3]); 
        dKdalpha2dalpha2 = dKdalpha2 .* (1 + permute(scaledDelta, [1 2 4 3])); 
        Khess = cat(4, cat(3, dKdalpha1dalpha1, dKdalpha1dalpha2), ...
            cat(3, dKdalpha2dalpha1, dKdalpha2dalpha2)); 
        
    end
end

end
