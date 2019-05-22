function [h]=hyp1f1_taylor2(a,b,z,tol)

if nargin < 4
    tol = 1e-6;
end

% Initialise a1, vector of individual terms, and b1, which stores the sum
% of the computed terms up to that point
[a,b,z] = expandSizes(3,a,b,z); 
iter = 0;
oldh = ones(size(z));
oldoldh = zeros(size(z));
h = ones(size(z)); 
while iter < 500
    r = (a+iter)./((iter+1)*(b+iter));
    h = oldh + (oldh - oldoldh).*r.*z; 
%     [r h]
    if all((abs((h - oldh)./h) < tol) & (abs((oldh - oldoldh)./h) < tol))
        break
    end
    oldoldh = oldh;
    oldh = h;
    iter = iter + 1; 
end

disp(iter)
