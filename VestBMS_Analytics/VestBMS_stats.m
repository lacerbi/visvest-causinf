function varargout = VestBMS_stats(n,mfits,bms)
%VESTBMS_STATS Compute statistics
%   VESTBMS_STATS(N) compute the N-th analysis.

% Luigi Acerbi 2016

% Load model fit results
if nargin < 2 || isempty(mfits)
    mfits = load('VestBMS_modelfits');
end
if nargin < 3; bms = []; end

mfits.mbag_uni.legend = 'Unisensory measurement';
mfits.mbag_biml.legend = 'Bisensory implicit CI';
mfits.mbag_bimu.legend = 'Bisensory explicit CI';
mfits.mbag_joint.legend = 'Joint datasets';

metric = 'loocv';
subjs = 1:11;

switch(n(1))
    
    case 1  % Compute goodness of fit
        % Usage: [R2,modelnames] = VestBMS_stats([1,k],mfits)
        
        switch(n(2))
            case 1
                fprintf('Explicit causal inference.\n');
                task = 3; bagidx = 2;
                flags = [0 1 0 0];
                mbag = mfits.mbag_bimu;
            case 2
                fprintf('Implicit causal inference.\n');
                task = 2; bagidx = 1;
                flags = [1 0 0 0];
                mbag = mfits.mbag_biml;
            case 3
                fprintf('Joint fits. Explicit causal inference.\n');
                task = [2 3 4]; bagidx = 3;
                flags = [0 0 1 0];
                mbag = mfits.mbag_joint;
                
        end
        
        if isempty(bms)
            [bms,mfits] = VestBMS_modelComparison(metric,1,mfits,[],[],flags);
        end
        modelnames = bms{bagidx}.modelnames;
        Nmodels = numel(modelnames);

        gof = [];
        modelsummary = ModelWork_summary(mbag);
        for i = 1:Nmodels
            mfit = ModelBag_get(mbag, modelsummary.dataid, readmodel(modelnames{i}, modelsummary), modelsummary.cnd);        
            if isempty(gof); gof = zeros(Nmodels,numel(mfit)); end
            [gof(i,:),gbar,H(i,:)] = VestBMS_gof(mfit,task);
            fprintf('Model %8s: pxp %.3f, R2 = %.3f ± %.3f.\n', modelnames{i}, bms{bagidx}.pxp(i), mean(gof(i,:)), stderr(gof(i,:)));
        end        
        fprintf('Model %8s: pxp %.3f, R2 = %.3f ± %.3f.\n', 'HIST', NaN, mean(gbar(1,:)), stderr(gbar(1,:)));
            
        varargout{1} = gof;
        varargout{2} = modelnames;
        varargout{3} = H;

    case 2 % Compute factor results
        
        flags = [1 1 1 1];
        if isempty(bms)
            [bms,mfits] = VestBMS_modelComparison(metric,1,mfits,[],[],flags);
            close(gcf);
        end
        
        % Results of explicit inference
        printbms(bms{2},'Explicit inference:');
        fprintf('%s is ~ %.2f times less likely to be the most representative model than the others.\n', bms{2}.factor{2}.factornames{3}, sum(bms{2}.factor{2}.tab(1:2))/bms{2}.factor{2}.tab(3));
        fprintf('%s is ~ %.2f times more likely than the others.\n', bms{2}.factor{2}.factornames{1}, sum(bms{2}.factor{2}.tab(1))/sum(bms{2}.factor{2}.tab(2:3)));
        
        printbms(bms{1},'Implicit inference:'); 
        fprintf('Ratio %s to %s ~ %.2f.\n', bms{1}.factor{1}.factornames{2}, bms{1}.factor{1}.factornames{1}, bms{1}.factor{1}.tab(2)/bms{1}.factor{1}.tab(1));
        fprintf('Ratio %s to %s ~ %.2f.\n', bms{1}.factor{3}.factornames{1}, bms{1}.factor{3}.factornames{2}, bms{1}.factor{3}.tab(1)/bms{1}.factor{3}.tab(2));
        
        printbms(bms{3},'Joint inference:'); 
        fprintf('Causal inference is ~ %.2f times more likely than forced fusion.\n', sum(bms{3}.factor{2}.tab(1:2))/sum(bms{3}.factor{2}.tab(3:4)));
        fprintf('Ratio %s to %s ~ %.2f.\n', bms{3}.factor{1}.factornames{2}, bms{3}.factor{1}.factornames{1}, bms{3}.factor{1}.tab(2)/bms{3}.factor{1}.tab(1));
        fprintf('Ratio %s to %s ~ %.2f.\n', bms{3}.factor{3}.factornames{1}, bms{3}.factor{3}.factornames{2}, bms{3}.factor{3}.tab(1)/bms{3}.factor{3}.tab(2));
        
        
    case 3  % Compute posteriors and overlap metrics
        % Usage: C = VestBMS_stats(3,mfits,bms)
        
        flags = [1 1 1 1];
        if isempty(bms)
            [bms,mfits] = VestBMS_modelComparison(metric,1,mfits,[],[],flags);
            close(gcf);
        end
        mfits = VestBMS_fixModels(mfits);
        
        fprintf('\nJoint fits (all datasets):\n');
        mbag_three = {mfits.mbag_biml, mfits.mbag_bimu, mfits.mbag_uni};
        modelsummary_three = {mfits.modelsummary_biml, mfits.modelsummary_bimu, mfits.modelsummary_uni};
        bms_three = {bms{1},bms{2},bms{4}};
        params = {'sigma_vest','sigma_vis_low','sigma_vis_med','sigma_vis_high','w_vest','w_vis','lambda'};
        for i = 1:numel(params)
            subplot(3,3,i);
            [postpdf,pdf,comp.C3(i,:),comp.P3(i,:,:),comp.pxp3(i,:,:)] = ModelPlot_parameterDistribution(mbag_three,params{i},bms_three,[],modelsummary_three,0);
            fprintf('%20s P(H0|data) = %.3f\n', [params{i} ','], comp.pxp3(i,1));
        end

        fprintf('\nJoint fits (bisensory datasets only):\n');
        mbag_two = {mfits.mbag_biml, mfits.mbag_bimu};
        bms_two = {bms{1},bms{2}};
        modelsummary_two = {mfits.modelsummary_biml, mfits.modelsummary_bimu};
        params = {'sigma_vest','sigma_vis_low','sigma_vis_med','sigma_vis_high','w_vest','w_vis','lambda','pcommon','kcommon','priorsigma','priorsigmadelta'};
        for i = 1:numel(params)
            subplot(3,4,i);
            [postpdf,pdf,comp.C2(i,:),comp.P2(i,:,:),comp.pxp2(i,:,:)] = ModelPlot_parameterDistribution(mbag_two,params{i},bms_two,[],modelsummary_two,0);
            fprintf('%20s P(H0|data) = %.3f\n', [params{i} ','], comp.pxp2(i,1));
        end
        
%         fprintf('\nJoint fits (unity judgements and unisensory only):\n');
%         mbag_two = {mfits.mbag_bimu, mfits.mbag_uni};
%         bms_two = {bms{2},bms{4}};
%         modelsummary_two = {mfits.modelsummary_bimu, mfits.modelsummary_uni};
%         params = {'sigma_vest','sigma_vis_low','sigma_vis_med','sigma_vis_high','w_vest','w_vis','lambda'};
%         for i = 1:numel(params)
%             [postpdf,pdf,comp.C2bis(i,:),comp.P2bis(i,:,:),comp.pxp2bis(i,:,:)] = ModelPlot_parameterDistribution(mbag_two,params{i},bms_two,[],modelsummary_two,0);
%             fprintf('%20s P(H0|data) = %.3f\n', [params{i} ','], comp.pxp2bis(i,1));
%         end
        
        varargout{1} = comp;
        
        
    case 4  % Compute rm-ANOVAs
        
        load('VestBMS_data.mat');
        
        switch(n(2))
            case 1  % Explicit inference (unity judgments)
                for i = subjs
                    for iNoise = 1:3
                        temp = data{i}.X.bimdisp{iNoise}{3};
                        [~,~,idx] = unique(temp(:,3));
                        for j = 1:max(idx)
                            X(i,j,iNoise) = sum(temp(idx == j,4) == 1) / sum(idx == j);
                        end                        
                    end
                end
                factors = {'noise', 'disparity'};
                tab = teg_repeated_measures_ANOVA([X(:, :, 1) X(:, :, 2) X(:, :, 3)], [3 size(X,2)], factors);
                rmANOVA_print(tab);        
                
            case 2  % Implicit inference (vestibular PSE)
                for i = subjs; X(i,:,:) = data{i}.psyright_mu'; end
                factors = {'noise', 's_vis'};
                tab = teg_repeated_measures_ANOVA([X(:, :, 1) X(:, :, 2) X(:, :, 3)], [3 size(X,2)], factors);
                rmANOVA_print(tab);        
                
            case 3  % Implicit inference (vestibular PSE)
                for i = subjs; X(i,:,:) = data{i}.psyright_mu'; end
                X(:,2,:) = -X(:,end,:);
                X = squeeze(nanmean(X(:,[1 2],:),2));
                [h,p,~,tab] = ttest(X)
                mean(X)./(tab.sd)
        end
        
        varargout{1} = tab;
        
end




end

%--------------------------------------------------------------------------
function model = readmodel(model, modelsummary)
%READMODEL Convert model in any form to model vector   
    if ischar(model)
        idx = find(strcmp(model,modelsummary.modelnames));
        if isempty(idx) || ~isscalar(idx)
            error(['Cannot find unique match for model ' model ' in model summary.']);
        end
        model = modelsummary.models(idx,:);
    elseif isscalar(model)
        model = modelsummary.models(model,:);        
    end
end

%--------------------------------------------------------------------------
function printbms(bms,txt)
%PRINTBMS Print results of Bayesian Model Selection

fprintf('\n%s\n\n',txt);

% Print factors
factor = bms.factor;
for i = 1:numel(factor)
    for j = 1:numel(factor{i}.tab)
        fprintf('%s %.2f\t', factor{i}.factornames{j},factor{i}.tab(j));
    end
    fprintf('\n');
end

% Print models
fprintf('Models (pxp):\n');
[y,ord] = sort(bms.pxp,'descend');
for i = 1:numel(y)
    fprintf('%s %.2f\n',bms.modelnames{ord(i)},y(i));
end
fprintf('\n');


end

%--------------------------------------------------------------------------
function rmANOVA_print(tab)

Nf = size(tab.R,1);
len = max(cellfun(@numel,tab.labels));
formatstr = ['%-' num2str(len) 's\tF(%.2f,%.2f) = %.2f,\teps = %.2g,\tp = %.3g,\teta^2_p = %.2f\n'];

for i = 1:Nf
    R = tab.R(i,:);
    eps0 = tab.eps0(i);
    fprintf(formatstr, tab.labels{i}, R(2)*eps0, R(3)*eps0, R(1), eps0, R(4), R(7)); 
end

end