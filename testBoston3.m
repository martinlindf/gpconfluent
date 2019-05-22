
dat = csvread('BostonHousing.csv',1,0);
load hermquad


cpdf = @(x,mu,sigma2,a,b) gammac(a+.5).*gammac(a+b)./(gammac(a).*gammac(a+b+.5).*sqrt(2*pi*sigma2)) .* hyp1f1g(a+.5,a+b+.5, -(x-mu).^2./(2*sigma2)); 
tpdf = @(x,mu,sigma2,nu) gammac((nu+1)/2)./(gammac(nu/2).*sqrt(pi*sigma2*nu)) .* (1+(x-mu).^2./(sigma2*nu)).^(-(nu+1)/2); 

desc = {'GP, Confluent', 'GP, Student''s t', 'GP, Gaussian'};

nTrain = 200; 
nFolds = 20; 

x = dat(:, 1:13);
x = (x - mean(x))./std(x);
y = dat(:, 14);  
y = (y - mean(y))./std(y);  
N = length(y);
rmseTest = zeros(nFolds, 2);
llTest = zeros(nFolds, 2); 
rmseFunTrain = zeros(nFolds, 2);
llFunTrain = zeros(nFolds, 2);
rmseFunTest = zeros(nFolds, 2);
llFunTest = zeros(nFolds, 2);


now = datestr(clock, 'yy-mm-dd_HH-MM-SS');
pool = gcp('nocreate');
if isempty(pool)
     pool = parpool([2 10]); 
end



parfor iFold = 1:nFolds
    idx = randperm(N);
    for iMethod = 1:2

        idxTrain = idx(1:nTrain); 
        xTrain = x(idxTrain,:);
        yTrain = y(idxTrain,:);


        idxTest = idx(nTrain+1:end); 
        xTest = x(idxTest, :);
        yTest = y(idxTest, :); 


        opt = struct; 
        switch iMethod
            case 1
                opt.noiseType = 'b'; 
            case 2
                opt.noiseType = 't'; 
        end

        opt.kerName = 'squaredExponentialnD'; 
        [par,hyp] = varGP(yTrain,xTrain,opt);
        


        Ktot = hyp.kernel.K([xTrain; xTest], hyp.alpha); 
        Ktot = Ktot+hyp.fudge*eye(size(Ktot)); 
        Kcross = Ktot(nTrain+1:end, 1:nTrain); 
        K = Ktot(1:nTrain, 1:nTrain);
        Kstar = Ktot(nTrain+1:end, nTrain+1:end); 




        MM = Kcross/((K));
        mpost = MM*yTrain;
        Kpost = Kstar - MM*Kcross';



        rmseTest(iFold,iMethod) = sqrt(mean((mpost - yTest).^2)); 

        switch iMethod
            case 1
                integrand = cpdf(mpost+sqrt(2*diag(Kpost)).*grid, yTest, hyp.sigma2, hyp.a, hyp.b);
                llTest(iFold,iMethod) = sum(log(1/sqrt(pi)*sum(weights.* integrand, 2))); 
            case 2
                integrand = tpdf(mpost+sqrt(2*diag(Kpost)).*grid, yTest, hyp.sigma2, hyp.nu); 
                llTest(iFold,iMethod) = sum(log(1/sqrt(pi)*sum(weights.* integrand, 2)));
        end
    end

end

delete(pool); 



save('-v7.3', sprintf('out/concreteResults-%s.mat', now), 'rmseTest', 'llTest','nTrain', 'nFolds', 'desc'); 

%%

