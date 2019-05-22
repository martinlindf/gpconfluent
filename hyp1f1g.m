function [h]=hyp1f1g(a,b,z,tol)

if nargin < 4
    tol = 1e-12;
end
[a,b,z] = expandSizes(3,a,b,z); 

asympLim = (abs(z)>150+1*abs(a)+1*abs(b));
case1 = (real(z)>=0) & ~asympLim; 
case2 = (real(z)<0) & ~asympLim;
case3 = asympLim; 
% pos = 1;
h = zeros(size(z));
h(case1) = hyp1f1_taylor(a(case1),b(case1),z(case1),tol);
h(case2) = exp(z(case2)).*hyp1f1_taylor(b(case2)-a(case2),b(case2),-z(case2), tol);
h(case3) = hyp1f1_asymp(a(case3),b(case3),z(case3),tol); 
