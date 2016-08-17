function [ctab,bias,se,mtab,stab] = VestBMS_showModelRecovery(task,mbag)
%VESTBMS_SHOWMODELRECOVERY Show model recovery results.

if nargin < 2; mbag = []; end
modelsummary = [];

% Datasets
Nsubjs = 11;
N = 10;
dataids = [];
for i = 1:4
    temp = bsxfun(@plus, (1:5)', ((1:Nsubjs)-1)*N) + (i-1)*(N*Nsubjs);
    dataids = [dataids; temp(:)];
end

switch lower(task(1))
    case 'u'; fitname = 'VestBMS_21001'; dataset = 'VestBMS_fakedata_unity';
    case 'l'; fitname = 'VestBMS_20001'; dataset = 'VestBMS_fakedata_bimloc';
    otherwise
        error('Unknown task. Specify either (U)nity judgment or bimodal (L)ocalization.');        
end
load(dataset);
if isempty(mbag)
    load(['.' filesep 'modelrecovery' filesep fitname]);
end
if isempty(modelsummary)
    modelsummary = ModelWork_summary(mbag);
end
[ctab,bias,se,mtab,stab] = ModelWork_plotModelRecovery(data(dataids),mbag,modelsummary);

end