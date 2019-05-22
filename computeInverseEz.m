function k = computeInverseEz(a, b, Ez)
% Solves Ez = (alpha./(alpha+beta)) .* hyp1f1g(alpha+1, alpha+beta+1, k) ./ hyp1f1g(alpha, alpha+beta, k) 
                % where alpha = hyper.a + 1/2, beta = hyper.b


assert(b==0.1, 'wrong value of b: only b = 0.1 is implemented! but it appears fairly insensitive.');

[agrid, Ezgrid, kgrid] = lookupTableInverseEz; 


k = interp2(agrid,Ezgrid,kgrid,a,Ez); 




