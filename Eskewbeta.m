
function [Ew, Elogw, Elog1mw, logNormConst] = Eskewbeta(alpha, beta, k)

lognormalizer = @(aa,bb,kk) log(gammac(aa)) + log(gammac(bb)) - log(gammac(aa+bb)) + log(hyp1f1g(aa,aa+bb, kk)); 


epsi = 1e-7; % error \approxeq epsi^2

Ew = (alpha./(alpha+beta)) .* hyp1f1g(alpha+1, alpha+beta+1, k) ./ hyp1f1g(alpha, alpha+beta, k) ; 
if nargout > 1
    Elogw = imag(lognormalizer(alpha+epsi*1i, beta, k)/epsi); 
    Elog1mw = imag(lognormalizer(alpha, beta+epsi*1i, k)/epsi); 
    logNormConst = lognormalizer(alpha,beta,k); 
end
 