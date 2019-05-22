function [param_k, hyper_k,y,x] = varGP(y,x,pars)
%% function [param_k, hyper_k,y,x] = varGP(y,x,pars)
% INPUT: y, vector of output variables
%        x, input of regressors
%        pars, struct of options, see below
% OUTPUT: 
% param_k, struct of variational parameters describing q_F(F) and q_Z(Z)
% hyper_k, struct of variational hyper-parameters (eta_Z)
%
%
% pars.iterMax / pars.iterMin : maximum/minimum # of variational iterations
% pars.iterProposal: number of initialization iterations
% pars.tol: ELBO optimization tolerance
% pars.noiseType: String setting noise type, default: 'b' for G-confluent, optional: 't' for Student's t
% par.kerName: kernel name, default: squaredExponentialnD

if nargin == 0
    % Default data
    rng(1)
    x = randn(100, 1)*.5+1.5*sign(rand(100,1)-0.5); x = sort(x);  
    y = trnd(1, 100, 1)*0.1+sin(2*x); 
    N = length(y);
end


% define default variational parameters
iterMax = 5000; 
iterMin = 100; 
iterVarMax = 1; 
iterProposal = 5; 
initThetaMu = [];
initThetaCov = [];
initSigma2 = []; 

tol = 1e-4; 

trueguess = 0; 

kerName = 'squaredExponentialnD'; 



noiseType = 'b';  

% noiseType = 't';

% noiseType = 'g'; 

N = length(y);


if nargin >= 3 && isstruct(pars)
    fn = fieldnames(pars); 
    for ii = 1:length(fn)
        eval([fn{ii} '= pars.' fn{ii} ';']); 
    end
end

%%







tmp=fitrgp(x,y, 'kernelfunction', 'squaredexponential'); 

% compute initialization guess using GP regression with Gauss Lik (matlab)


thRob = resubPredict(tmp);
sigma2Hat = tmp.Sigma.^2*10; 
cholprecHat = eye(N); 

if trueguess
    thRob = 5*sin(2*pi*x);
end

switch lower(kerName)
case 'squaredexponential'
    ker = struct('name', 'squaredExponential', ...
                'alphanames', ['logscale '; 'loginvstd'], ...
                'alphaInit', combvec(-3:3:3, -2:2:2)', ...
                'alphaLB', [-5 -10], ...
                'alphaUB', [+5 +10], ...
                'alphaCon', [], ...
                'alphaTarget', [],... % left empty -> compute automatically. Otherwise, include gradient and Hessian. 
                'K', @kernelSE ...
                );
case 'squaredexponentialnd'
    nx = size(x,2); 
    ker = struct('name', 'squaredExponentialnD', ...
                'alphanames', ['logscale '; 'loginvstd'], ...
                'alphaInit', [(-3:3:3)', zeros(3,nx)], ...
                'alphaLB', [-5 -10*ones(1,nx)], ...
                'alphaUB', [+5 +10*ones(1,nx)], ...
                'alphaCon', [], ...
                'alphaTarget', [],... % left empty -> compute automatically. Otherwise, include gradient and Hessian. 
                'K', @kernelSEnD ...
                );

case 'custom'
    ker = customKerStruct; 
end



hyper_0 = struct(...'kernel', ker,...
                'sigma2', sigma2Hat, ...
                'nu', 5, ...
                'a', 5/2, ...
                'b', 0.1, ...
                'alpha', ker.alphaInit(1,:)', ... 
                'fudge', 1e-8, ...
                'type', noiseType ...
                  ); 
              
              
hyper_0.kernel = ker; 



param_0 = struct('f', ...
                    struct('mu', thRob, ...
                        'cholprec', cholprecHat,...
                        'cholcov', inv(cholprecHat)), ...
                'z', ...
                    struct('a', hyper_0.a, 'b', hyper_0.b, 'k', 0)); 


param_0.f.cov = param_0.f.cholcov * param_0.f.cholcov'; 

tic
nu_grid = [2; 4; 6; 8];
sigma2_grid = sigma2Hat*[ .01 .03 .1 .3 1 3 10 30 100]; 
% 
elboInit = -inf(length(nu_grid), size(hyper_0.kernel.alphaInit,1), length(sigma2_grid)); 
for ii = 1:length(nu_grid)
    for jj = 1:size(hyper_0.kernel.alphaInit,1)
        for kk = 1:length(sigma2_grid)
        
        
            hypInit{ii,jj,kk} = hyper_0;  %#ok<AGROW>
            hypInit{ii,jj,kk}.nu = nu_grid(ii); %#ok<AGROW>
            hypInit{ii,jj,kk}.a = (nu_grid(ii))/2; %#ok<AGROW>
            hypInit{ii,jj,kk}.sigma2 = sigma2_grid(kk);  %#ok<AGROW>
            hypInit{ii,jj,kk}.alpha = hypInit{ii,jj,kk}.kernel.alphaInit(jj,:)'; %#ok<AGROW>
        
            EZinit = 0.9*ones(N,1);  
            
            if hypInit{ii,jj,kk}.type == 't'
                EZinit = (EZinit / mean(EZinit)); % typically EZ = 1 on average
                % Solve EZ = a/b where a = (hyper.nu + 1)/2 -> b = a/EZ
                Zinit = struct('a', (hypInit{ii,jj,kk}.nu + 1)/2, 'b', ((hypInit{ii,jj,kk}.nu + 1)/2)./EZinit, 'k', 0); 
            elseif hypInit{ii,jj,kk}.type == 'b' 
                EZinit = (EZinit / mean(EZinit)) * hypInit{ii,jj,kk}.a / (hypInit{ii,jj,kk}.a + hypInit{ii,jj,kk}.b); 
                % if possible, adjust so EZ is a/(a+b) on average
                EZinit = (EZinit / max(EZinit))*0.9; % but EZ < 1 so this is most important. 
                EZinit = max(EZinit, 0.01);  % lookup table implemented in 0.01 ... 0.9.
                
                
                kinit = computeInverseEz(hypInit{ii,jj,kk}.a, hypInit{ii,jj,kk}.b, EZinit); 
                
                
                
                Zinit = struct('a', hypInit{ii,jj,kk}.a + 1/2, 'b', hypInit{ii,jj,kk}.b, 'k', kinit); 
                
                
            elseif hypInit{ii,jj,kk}.type == 'g'
                Zinit = struct('a', 1,'b',ones(N,1)); 
            else 
                error('unknown noise model!'); 
            end
            
            param_0 = struct('f', ...
                                struct('mu', thRob, ...
                                    'cholprec', cholprecHat,...
                                    'cholcov', inv(cholprecHat)), ...
                            'z', ...
                                Zinit); 
            param_0.f.cov = param_0.f.cholcov * param_0.f.cholcov'; 


            param_0 = updateParamStatsF(param_0); 
            param_0 = updateParamStatsZ(param_0, hyper_0.type); 
            
            try
                [elboInit(ii,jj,kk), paramTunedInit{ii,jj,kk}, hyperTunedInit{ii,jj,kk}] = ...
                    run_var_GP(y, x, param_0, hypInit{ii,jj,kk}, iterProposal, iterProposal, 1, tol); %#ok<AGROW>
            catch
                elboInit(ii,jj,kk) = -inf; 
                paramTunedInit{ii,jj,kk} = param_0;  %#ok<AGROW>
                hyperTunedInit{ii,jj,kk} = hypInit{ii,jj,kk};  %#ok<AGROW>
            end
        end

    end
end

[~, init_ind] = max(elboInit(:));
param_0 = paramTunedInit{init_ind}; 
hyper_0 = hyperTunedInit{init_ind}; 

param_0 = updateParamStatsF(param_0); 
param_0 = updateParamStatsZ(param_0, hyper_0.type); 


            
            
            
[elbo, param_k, hyper_k, iterMax, elbo_hist, param_hist, hyper_hist] = ...
    run_var_GP(y, x, param_0, hyper_0, iterMin, iterMax, iterVarMax, tol);  %#ok<ASGLU>



if(toc>1)
    toc
end

%%


% 
if hyper_k.type == 't'
    str = '                    sigma2         nu ';
    for ii = 1:size(ker.alphanames,1)
        str = [str '      ' ker.alphanames(ii,:)]; %#ok<AGROW>
    end

    disp(str); 
    fprintf('Initial hyper:')
    disp([hyper_0.sigma2, hyper_0.nu, hyper_0.alpha(:)']);
    fprintf('Final hyper:  '); 
    disp([hyper_k.sigma2, hyper_k.nu, hyper_k.alpha(:)']);
    fprintf('\n'); 
elseif hyper_k.type == 'b'
    str = '                    sigma2         a           b  ';
    for ii = 1:size(ker.alphanames,1)
        str = [str '      ' ker.alphanames(ii,:)]; %#ok<AGROW>
    end

    disp(str); 
    fprintf('Initial hyper:')
    disp([hyper_0.sigma2, hyper_0.a, hyper_0.b, hyper_0.alpha(:)']);
    fprintf('Final hyper:  '); 
    disp([hyper_k.sigma2, hyper_k.a, hyper_k.b, hyper_k.alpha(:)']);
    fprintf('\n'); 
elseif hyper_k.type == 'g'
    str = '                    sigma2         ';
    for ii = 1:size(ker.alphanames,1)
        str = [str '      ' ker.alphanames(ii,:)]; %#ok<AGROW>
    end

    disp(str); 
    fprintf('Initial hyper:')
    disp([hyper_0.sigma2, hyper_0.alpha(:)']);
    fprintf('Final hyper:  '); 
    disp([hyper_k.sigma2, hyper_k.alpha(:)']);
    fprintf('\n'); 
else
    warning('unknown noise type!'); 
end


end

