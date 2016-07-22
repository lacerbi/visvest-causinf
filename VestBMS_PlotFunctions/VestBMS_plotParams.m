function VestBMS_plotParams(nid)

fitnames = {'uni','biml','bimu','joint'};
if 0
    modeln = [2 1 1 1];
else
    modeln = [2 1 3 1];    
end

mbags = load('VestBMS_modelfits.mat');

ids = 1:11; % Human only

for i = 1:numel(fitnames)
    mbag = mbags.(['mbag_' fitnames{i}]);
    modelsummary = mbags.(['modelsummary_' fitnames{i}]);
    mfits{i} = ModelBag_get(mbag,modelsummary.dataid(ids,:),modelsummary.models(modeln(i),:),modelsummary.cnd);
end

figure;
ModelPlot_plotParameters('VestBMS',...
    mfits{1}{nid},mfits{2}{nid},mfits{3}{nid},mfits{4}{nid})


end