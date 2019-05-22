
function bound = ELBO_GP(param, hyper, y, x) 


[N,n] = size(x);
    


ElogpriorF = - 0.5 * alphaTargetGeneric(hyper, param.stats.EfSq, x, hyper.alpha, n) - n/2*log(2*pi);
% z prior
if hyper.type == 't'
    ElogpriorZ = (hyper.nu/2-1) * param.stats.EsumlogZ - hyper.nu/2 * param.stats.EsumZ ...
                + N * (hyper.nu/2 * log(hyper.nu/2) - gammaln(hyper.nu/2) );
elseif hyper.type == 'b'
    ElogpriorZ = (hyper.a-1)*param.stats.EsumlogZ + (hyper.b-1)*param.stats.Esumlog1mZ - N*betaln(hyper.a, hyper.b); 
elseif hyper.type == 'g'
    ElogpriorZ = 0; 
else 
    error('unknown noise model!'); 
end
% likelihood
Eloglik = - 1/(2*(hyper.sigma2)) * (param.stats.EZ' * ((y-param.f.mu).^2 + diag(param.f.cov))/N) ...
     - N/2 * log(2*pi) + 0.5*param.stats.EsumlogZ - N/2 * log(hyper.sigma2);
 
% we have computed all components, so the ELBO is...
bound = ElogpriorF + ElogpriorZ + Eloglik...
                + param.stats.entropyF + param.stats.entropyZ; 
            
end
