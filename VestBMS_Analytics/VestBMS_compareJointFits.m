%VESTBMS_COMPAREJOINTFITS Compare joint vs separate fits.
function [msep,mjoint] = VestBMS_compareJointFits(type)

mfits = load('VestBMS_modelFits.mat');
subjs = 1:11;     % Humans only

switch lower(type(1))
    case 'f'; modelnames = {'BP','CXD','CX','CXD'};
    case 'b'; modelnames = {'BP','BPD','BPD','BPD'};
    otherwise; error('Available type are (B)ayesian and (F)ixed criterion.');
end

if 0
    tasks = {'uni','biml','bimu','joint'};
    display('Comparing JOINT fits to unisensory + bimodal localization + unity judgments.');
else
    modelnames(1) = [];
    tasks = {'biml','bimu','semi'};
    display('Comparing SEMI-JOINT fits to bimodal localization + unity judgments.');
end

metrics = {'aic','aicc','bic','dic','loocv','maploglike','marginallike_rlr','marginallike_whmg','marginallike_whmu'};

for j = 1:numel(metrics)

    score = NaN(numel(subjs),numel(tasks));
    for i = 1:numel(tasks)
        m = mfits.(['modelsummary_',tasks{i}]);
        idx = find(strcmp(m.modelnames,modelnames{i}));
        if isempty(idx) || ~isscalar(idx)
            error(['Cannot find a single model for model name ' modelnames{i} ' in ' ['modelsummary_',tasks{i}] '.']);
        end
        if isfield(m,metrics{j}) && ~isempty(m.(metrics{j}))
            score(:,i) = m.(metrics{j})(subjs,idx);
        end
    end

    msep.(metrics{j}) = sum(score(:,1:end-1),2);
    mjoint.(metrics{j}) = score(:,end);
end

[msep,mjoint]

end
