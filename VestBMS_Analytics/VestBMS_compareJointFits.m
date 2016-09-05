%VESTBMS_COMPAREJOINTFITS Compare joint vs separate fits.
function [msep,mjoint] = VestBMS_compareJointFits(type)

mfits = load('VestBMS_modelFits.mat');
subjs = 1:11;     % Humans only

switch lower(type(1))
    case 'f'; modelnames = {'BP','CXD','CX','CXD'};
    case 'b'; modelnames = {'BP','BPD','BPD','BPD'};
    otherwise; error('Available type are (B)ayesian and (F)ixed criterion.');
end
tasks = {'uni','biml','bimu','joint'};

metrics = {'aic','aicc','bic','dic','loocv','maploglike','marginallike_rlr','marginallike_whmg','marginallike_whmu'};

for j = 1:numel(metrics)

    score = NaN(numel(subjs),numel(tasks));
    for i = 1:numel(tasks)
        m = mfits.(['modelsummary_',tasks{i}]);
        idx = find(strcmp(m.modelnames,modelnames{i}));
        if isempty(idx) || ~isscalar(idx)
            error(['Cannot find a single model for model name ' modelnames{i} ' in ' ['modelsummary_',tasks{i}] '.']);
        end
        score(:,i) = m.(metrics{j})(subjs,idx);
    end

    msep.(metrics{j}) = sum(score(:,1:3),2);
    mjoint.(metrics{j}) = score(:,4);
end

[msep,mjoint]

end
