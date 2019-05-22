
now = datestr(clock, 'yy-mm-dd_HH-MM-SS');
pool = gcp('nocreate');
if isempty(pool)
     pool = parpool([2 12]); 
end
 
load hermquad

cpdf = @(x,mu,sigma2,a,b) gammac(a+.5).*gammac(a+b)./(gammac(a).*gammac(a+b+.5).*sqrt(2*pi*sigma2)) .* hyp1f1g(a+.5,a+b+.5, -(x-mu).^2./(2*sigma2)); 
tpdf = @(x,mu,sigma2,nu) gammac((nu+1)/2)./(gammac(nu/2).*sqrt(pi*sigma2*nu)) .* (1+(x-mu).^2./(sigma2*nu)).^(-(nu+1)/2); 
npdf = @(x,mu,sigma2) 1./(sqrt(2*pi*sigma2)) .* exp(-0.5*(x-mu).^2./sigma2);

numIter = 100; 
Ntrain = 100;
Ntest = 100;
N = Ntrain + Ntest; 
n = 10; 
numOut = 20; 

outScale = [1 3 10 30]; 
desc = {'GP, Confluent', 'GP, Student''s t', 'GP, Gaussian'};

numOutScale = length(outScale);
rmseFunTrain = zeros(2,numOutScale,numIter);
llFunTrain = zeros(2,numOutScale,numIter);
rmseFunTest = zeros(2,numOutScale,numIter);
llFunTest = zeros(2,numOutScale,numIter);
rmseTest = zeros(2,numOutScale,numIter);
llTest = zeros(2,numOutScale,numIter);
dat = repmat(struct('y',zeros(N,1), 'x', zeros(N,n), 'f0',zeros(N,1)), [numOutScale numIter]);

rng(1); 
parfor ii = 1:numIter
    for iOut = 1:numOutScale
        % generate data 
        x = rand(N,n); % uniformly on hypercube
        f0 = sin(pi*x(:,1).*x(:,2)) + 20*(x(:, 3)-0.5).^2 + 10*x(:,4) + 5*x(:,5);  % friedman's function 

        y0 = f0 + randn(N,1);
        out = randi(N, numOut, 1);
        y = y0;
        y(out) = 15+outScale(iOut)*randn(1,numOut); 

        dat(iOut, ii).y = y;
        dat(iOut, ii).x = x;
        dat(iOut, ii).f0 = f0;
        xTrain = x(1:Ntrain, :);
        yTrain = y(1:Ntrain,:); 
        f0Train = f0(1:Ntrain,:); 
        xTest = x(Ntrain+1:end, :);
        yTest = y(Ntrain+1:end,:);
        f0Test = f0(Ntrain+1:end, :);
        
        for iMethod = 1:2
            opt = struct; 
            switch iMethod
                case 1
                    opt.noiseType = 'b'; 
                case 2
                    opt.noiseType = 't';  
            end
            opt.kerName = 'squaredExponentialnD'; 


            

            [par,hyp] = varGP(yTrain,xTrain,opt);
            rmseFunTrain(iMethod, iOut, ii) = sqrt(mean((f0Train - par.f.mu).^2)); 
            llFunTrain(iMethod, iOut, ii) = -0.5*((f0Train - par.f.mu)'/par.f.cov)*(f0Train-par.f.mu) - log(det(par.f.cov))/2; 
            
            
            Ktot = hyp.kernel.K([xTrain; xTest], hyp.alpha); 
            Ktot = Ktot+hyp.fudge*eye(size(Ktot)); 
            Kcross = Ktot(Ntrain+1:end, 1:Ntrain); 
            K = Ktot(1:Ntrain, 1:Ntrain);
            Kstar = Ktot(Ntrain+1:end, Ntrain+1:end); 




            MM = Kcross/((K));
            mpost = MM*yTrain;
            Kpost = Kstar - MM*Kcross';

            rmseFunTest(iMethod, iOut, ii) = sqrt(mean((mpost - f0Test).^2));
            llFunTest(iMethod, iOut, ii) = -0.5*((f0Test - mpost)'/Kpost)*(f0Test - mpost)- log(det(Kpost))/2; 
            
            
            rmseTest(iMethod, iOut, ii) = sqrt(mean((mpost - yTest).^2)); 

            switch iMethod
                case 1
                    integrand = cpdf(mpost+sqrt(2*diag(Kpost)).*grid, yTest, hyp.sigma2, hyp.a, hyp.b);
                    llTest(iMethod, iOut, ii) = sum(log(1/sqrt(pi)*sum(weights.* integrand, 2))); 
                case 2
                    integrand = tpdf(mpost+sqrt(2*diag(Kpost)).*grid, yTest, hyp.sigma2, hyp.nu); 
                    llTest(iMethod, iOut, ii) = sum(log(1/sqrt(pi)*sum(weights.* integrand, 2)));
            end
        
        end
    end
end

delete(pool); 

% save tmpresults-backup
save('-v7.3', sprintf('out/friedmanResults-%s.mat', now), 'outScale', 'rmseFunTest', 'llFunTest', 'rmseFunTrain', 'llFunTrain', 'rmseTest', 'llTest', 'N', 'Ntrain', 'Ntest', 'n', 'numOut', 'dat', 'desc'); 

%%
%%

