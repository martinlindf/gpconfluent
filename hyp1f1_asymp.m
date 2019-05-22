function [h3]=hyp1f1_asymp(a,b,z,tol)

if nargin < 4
    tol = 1e-6;
end
% z = -z; 
% Initialise a1, vector of individual terms, and b1, which stores the sum
% of the computed terms up to that point
[a,b,z] = expandSizes(3,a,b,z); 
iter = 0;
term = ones(size(z)); 
h1 = ones(size(z)); 
iter = 0;
while iter <= 12
    term = term .* (b-a+iter) .* (1-a+iter) ./(1+iter) ./ z;
    h1 = h1 + term; 
%     [r h]
%     if all((abs(term./h1) < tol))
%         break
%     end
    iter = iter + 1; 
end
% iter
term = ones(size(z)); 
h2 = ones(size(z)); 
iter = 0;
while iter <= 12
    term = term .* (a+iter) .* (1+a-b+iter) ./(1+iter) ./ (-z);
    h2 = h2 + term; 
%     [r h]
%     if all((abs(term./h2) < tol))
%         break
%     end
    iter = iter + 1; 
end
% iter
phi = angle(z); 
secondsumSign = (phi > -0.5*pi & phi < 1.5*pi)*2-1; 

firstsum = exp(z).*z.^(a-b)./gammac(a).*h1; 
secondsum = exp(secondsumSign.*pi.*1i.*a).*z.^(-a)./gammac(b-a).*h2; 
% secondsumBranch2 = exp(+pi*1i*a).*z.^(-a)./gammac(b-a).*h2;

h3 = zeros(size(z)); 

allreal = (imag(a)==0 & imag(b)==0 & imag(z)==0); 
case1 = (allreal & real(z)>0);
case2 = (allreal & real(z)<0);
case3 = (~allreal & imag(z)==0);
case4 = (~allreal & imag(z)~=0);

h3(case1) = real(firstsum(case1));
h3(case2) = real(firstsum(case2) + secondsum(case2));
h3(case3) = firstsum(case3)+secondsum(case3);
h3(case4) = firstsum(case4) + secondsum(case4);

h3 = h3.*gammac(b); 


% if (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)>0)
%     h3 = real(exp(z)*z^(a-b)/gammac(a)*h1);
% elseif (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)<0)
%     h3 = real((exp(z)*z^(a-b)/gammac(a)*h1...
%         +exp(-pi*1i*a)*z^(-a)/gammac(b-a)*h2));
% elseif (imag(z)==0)
%     h3 = exp(z)*z^(a-b)/gammac(a)*h1;
% elseif (real(z)>0)
%     h3 = exp(z)*z^(a-b)/gammac(a)*h1...
%         +exp(pi*1i*a)*z^(-a)/gammac(b-a)*h2;
% elseif (real(z)<0)
%     h3 = exp(z)*z^(a-b)/gamfun(a)*h1...
%         +exp(-pi*1i*a)*z^(-a)/gammac(b-a)*h2;
% end
% 
% % Special case where a, b and z are real
% if imag(a)==0 && imag(b)==0 && imag(z)==0
%     h = real(h3);
% else
%     h = h3;
% end


% function [h] = hypfun_M_asymptotica(a,b,z,tol)
% % function [h] = hypfun_M_asymptotica(a,b,z,tol)
% % 
% % Computes the hypergeometric function \mathbf{M}(a;b;z), using an
% % asymptotic series (Method (a)), up to a tolerance tol.
% % 
% % Copyright John W. Pearson 2014
% 
% 
% % Sum series
% a1 = 1; b1 = 1;
% for j = 1:500
%     % Update a1(j) and b1 in terms of previously computed a1(j-1) and b1
%     a1 = (b-a+j-1)*(-a+j)/j/z*a1;
%     b1 = b1+a1;
%     % If stopping criterion is satisfied
%     if abs(a1)/abs(b1)<tol
%         break
%     end
%     % No convergence
%     if (j==500)
%         [' ' num2str(j) ' terms computed']
%         return
%     end
% end
% 
% % Sum series
% c1 = 1; d1 = 1;
% for k = 1:500
%     % Update c1 and d1
%     c1 = (a+k-1)*(a-b+k)/k/(-z)*c1;
%     d1 = d1+c1;
%     % Stopping criterion
%     if abs(c1)/abs(d1)<tol
%         break
%     end
%     % Specify if 500 terms computed
%     if (k==500)
%         [' ' num2str(k) ' terms computed']
%         return
%     end
% end
% 
% % Take last terms computed
% h1 = b1; h2 = d1;
% 
% % Compute an asymptotic series; which one depends on location of variables
% % in the complex plane
% if (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)>0)
%     h3 = real(exp(z)*z^(a-b)/gamfun(a)*h1);
% elseif (imag(a)==0 && imag(b)==0 && imag(z)==0 && real(z)<0)
%     h3 = real((exp(z)*z^(a-b)/gamfun(a)*h1...
%         +exp(-pi*1i*a)*z^(-a)/gamfun(b-a)*h2));
% elseif (imag(z)==0)
%     h3 = exp(z)*z^(a-b)/gamfun(a)*h1;
% elseif (real(z)>0)
%     h3 = exp(z)*z^(a-b)/gamfun(a)*h1...
%         +exp(pi*1i*a)*z^(-a)/gamfun(b-a)*h2;
% elseif (real(z)<0)
%     h3 = exp(z)*z^(a-b)/gamfun(a)*h1...
%         +exp(-pi*1i*a)*z^(-a)/gamfun(b-a)*h2;
% end
% 
% % Special case where a, b and z are real
% if imag(a)==0 && imag(b)==0 && imag(z)==0
%     h = real(h3);
% else
%     h = h3;
% end
