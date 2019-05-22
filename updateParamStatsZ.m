function param = updateParamStatsZ(param, type)
N = length(param.z.b); 
if type == 't'
    %update the (oversufficient) statistics.
    param.stats.EZ = param.z.a ./ param.z.b; % [N|1]
    param.stats.EsumZ = sum(param.stats.EZ); %[1|1]
    param.stats.EsumlogZ = sum( psi(param.z.a) - log(param.z.b)); % [1|1]
    param.stats.entropyZ = N*param.z.a - sum(log(param.z.b)) ...
        + N*gammaln(param.z.a) + N*(1-param.z.a)*psi(param.z.a);  % E -log q(z) [1|1]
elseif type == 'b' 
    [EZ, ElogZ, Elog1mZ, logNormConst] = Eskewbeta(param.z.a, param.z.b, param.z.k); 
    param.stats.EZ = EZ;
    param.stats.EsumZ = sum(param.stats.EZ); %[1|1]
    param.stats.EsumlogZ = sum(ElogZ); % [1|1]
    param.stats.Esumlog1mZ = sum(Elog1mZ); % [1|1]
%     [param.stats.EsumlogZ/length(ElogZ), param.stats.Esumlog1mZ/length(ElogZ), (param.stats.EsumlogZ + param.stats.Esumlog1mZ)/length(ElogZ)]
    param.stats.entropyZ = -sum((param.z.a-1) .* ElogZ + (param.z.b-1) .* Elog1mZ + param.z.k .* EZ - logNormConst);  % E -log q(z) [1|1]
elseif type == 'g'
    param.stats.EZ = ones(N,1);
    param.stats.EsumZ = N; 
    param.stats.EsumlogZ = 0; % [1|1]
    param.stats.entropyZ = 0; 
else 
    error('unknown noise model!'); 
end

% assert(