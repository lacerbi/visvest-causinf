function [R2,Y,Ypred] = VestBMS_gof(mfit,task)
%VESTBMS_GOF Compute model's goodness-of-fit.

type = 'entropy';
Ngen = 50;          % How many fake datasets for standard R2
data = [];
if ~iscell(mfit); mfit = {mfit}; end

R2 = zeros(1,numel(mfit));

for i = 1:numel(mfit)
    m = mfit{i};
    if isempty(m.mp)
        if isempty(data); load('VestBMS_data.mat'); end
        m = ModelWork_loadFields('VestBMS',m,data); 
    end
    
    switch lower(type)
        
        case 'generalized'        
            ll_0 = 0;
            ntot = 0;        

            for iTask = 1:numel(task)

                switch task(iTask)
                    case {1,2,4}
                        a = 1; b = -1; 
                    case 3
                        a = 1; b = 2;
                end

                switch task(iTask)
                    case {2,3}
                        X = []; Y = [];
                        for iNoise = 1:3; X = [X; m.X.bimodal{iNoise}{task(iTask)}]; end                
                    case 4
                        X = []; Y = [];
                        for iNoise = 1:4; X = [X; m.X.unimodal{iNoise}]; end
                end

                k = sum(X(:,end) == a);
                n = sum(X(:,end) == a | X(:,end) == b);
                ntot = ntot + n;
                p = k/n;
                ll_0 = ll_0 + k*log(p) + (n-k)*log(1-p);
            end

            ll_pred = m.metrics.maploglike;      
            R2(i) = 1 - exp(2/ntot*(ll_0-ll_pred));
            R2(i) = -(ll_pred-ll_0)/ll_0;
        
        case 'r2'
            m.sampling = []; % Ignore samples, use MLE
            gendata = VestBMS_gendata(Ngen,m);
            fakedata = VestBMS_analytics(gendata);


            for iTask = 1:numel(task)

                switch task(iTask)
                    case {1,2,4}
                        a = 1; b = -1; 
                    case 3
                        a = 1; b = 2;
                end

                Y = []; Ypred = [];
                for iNoise = 1:3
                    X = m.X.bimodal{iNoise}{3};
                    bins{iNoise} = unique(X(:,3) + 1000*X(:,4));
                    Y = [Y; computepred(X,bins{iNoise},a,b)];
                    G = [];
                    for ii = 1:Ngen
                        G = [G; fakedata{ii}.X.bimodal{iNoise}{3}];            
                    end
                    Ypred = [Ypred; computepred(G,bins{iNoise},a,b)];
                end
            end
            R2(i) = rsquare(Y,Ypred);
            
        case 'entropy'            
            x = []; D = [];
            % Unimodal data
            if any(task == 4)
                for iNoise = 1:4
                    D = [D; m.X.unibins{iNoise}];
                end
                x{end+1} = D;
            end
            % Bimodal data
            for iTask = task(task < 4)
                D = [];
                for iNoise = 1:3
                    D = [D; m.X.bimbins{iNoise}{iTask}];
                end
                x{end+1} = D;
            end
                        
            [R2(i),Y(i),output] = gofit(x,m.metrics.loocv,'gra');
            Ypred(i) = output.H / nansum(x{1}(:));
            % [output.H sqrt(output.Hvar)]
            % output
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