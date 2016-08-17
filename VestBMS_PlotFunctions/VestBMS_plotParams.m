function minq = VestBMS_plotParams(nids,modelfamily,joint,noplot)
%VESTBMS_PLOTPARAMS Plot and compare parameters from different tasks.

if nargin < 2 || isempty(modelfamily); modelfamily = 'bayes'; end
if nargin < 3 || isempty(joint); joint = 0; end
if nargin < 4 || isempty(noplot); noplot = 0; end

if joint
    fitnames = {'uni','biml','bimu','joint'};
else
    fitnames = {'uni','biml','bimu'};
end

switch(lower(modelfamily(1)))
    case 'b';           % Bayesian
        modeln = [2 10 6 2]; family = 'Bayesian';
    case {'f','c'};     % Fixed criterion
        modeln = [2 8 5 1]; family = 'Fixed';
    otherwise
        error('Unknown model family. Known families are (B)ayesian and (F)ixed criterion.');
end

fprintf('Compare parameters from the %s model family.\n', upper(family));

mbags = load('VestBMS_modelfits.mat');

ids = 1:11; % Human only

for i = 1:numel(fitnames)
    mbag = mbags.(['mbag_' fitnames{i}]);
    modelsummary = mbags.(['modelsummary_' fitnames{i}]);
    mfits{i} = ModelBag_get(mbag,modelsummary.dataid(ids,:),modelsummary.models(modeln(i),:),modelsummary.cnd);
end

for i = 1:numel(nids)
    if ~noplot; figure; end
    nid = nids(i);
    if numel(mfits) == 4
        minq(i,:) = ModelPlot_plotParameters('VestBMS',...
            mfits{1}{nid},mfits{2}{nid},mfits{3}{nid},mfits{4}{nid},noplot);
    else
        minq(i,:) = ModelPlot_plotParameters('VestBMS',...
            mfits{1}{nid},mfits{2}{nid},mfits{3}{nid},noplot);
    end
end


end