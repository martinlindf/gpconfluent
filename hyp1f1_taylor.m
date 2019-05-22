function [h]=hyp1f1_taylor(a,b,z,tol)

if nargin < 4
    tol = 1e-6;
end

% Initialise a1, vector of individual terms, and b1, which stores the sum
% of the computed terms up to that point
[a,b,z] = expandSizes(3,a,b,z); 
iter = 0;
oldterm = ones(size(z));
h = oldterm; 
while iter < 500
    term = oldterm.*(a+iter)./(b+iter).*z/(iter+1); 
    h = h + term; 
    if all((abs(term./h) < tol) & (abs(oldterm./h) < tol))
        break
    end
    oldterm = term; 
    iter = iter + 1; 
end

% disp(iter)
