function [R2,Y,Ypred] = VestBMS_gof(mfit,task)
%VESTBMS_GOF Compute model's goodness-of-fit.

generalized = 1;    % Compute generalized R2
Ngen = 50;          % How many fake datasets for standard R2
data = [];
if ~iscell(mfit); mfit = {mfit}; end

if task == 2; a = 1; b = -1; elseif task == 3; a = 1; b = 2; end

R2 = zeros(1,numel(mfit));

for i = 1:numel(mfit)
    m = mfit{i};
    if isempty(m.mp)
        if isempty(data); load('VestBMS_data.mat'); end
        m = ModelWork_loadFields('VestBMS',m,data); 
    end
    
    if generalized
        ll_pred = m.metrics.maploglike;        
        X = []; Y = [];
        for iNoise = 1:3; X = [X; m.X.bimodal{iNoise}{task}]; end
        k = sum(X(:,end) == a);
        n = sum(X(:,end) == a | X(:,end) == b);
        p = k/n;
        ll_0 = -binolike(p,n,p);

        R2(i) = 1 - exp(2/n*(ll_0-ll_pred));
    else
        m.sampling = []; % Ignore samples, use MLE
        gendata = VestBMS_gendata(Ngen,m);
        fakedata = VestBMS_analytics(gendata);
        
        Y = []; Ypred = [];
        for iNoise = 1:3
            X = m.X.bimodal{iNoise}{3};
            bins{iNoise} = unique(X(:,3) + 1000*X(:,4));
            Y = [Y; computepred(X,bins{iNoise})];
            G = [];
            for ii = 1:Ngen
                G = [G; fakedata{ii}.X.bimodal{iNoise}{3}];            
            end
            Ypred = [Ypred; computepred(G,bins{iNoise})];
        end
        R2(i) = rsquare(Y,Ypred);
    end
    

end



end

%--------------------------------------------------------------------------
function Y = computepred(X,bins,a,b)
    nb = numel(bins);
    xx = X(:,3) + 1000*X(:,4);
    Y = NaN(nb,1);
    for j = 1:nb
        yy = X(xx == bins(j),end);
        Y(j) = sum(yy == a)/sum(yy == a | yy == b);
    end
end