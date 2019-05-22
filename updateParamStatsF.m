function param = updateParamStatsF(param)
%update the (oversufficient) statistics.

cholcovAndMu = [param.f.cholcov param.f.mu]; 
n = length(param.f.mu); 

param.stats.EfSq = cholcovAndMu*cholcovAndMu'; % E theta theta' [n|n]
assert(all(isreal(param.stats.EfSq(:)))); 
param.stats.entropyF = n/2 * (1+log(2*pi)) - sum(log(abs(diag(param.f.cholprec))));  % E log q(f) [1|1]