%VESTBMS_COMPAREJOINTFITS Compare joint vs separate fits.
function [msep,mjoint] = VestBMS_compareJointFits()

[~,modelsummary{1}] = ModelWork_collectFits('VestBMS','uni*',[],[]);
[~,modelsummary{2}] = ModelWork_collectFits('VestBMS','bim-l*',[],[]);
[~,modelsummary{3}] = ModelWork_collectFits('VestBMS','bim-u*',[],[]);
[~,modelsummary{4}] = ModelWork_collectFits('VestBMS','joint*',[],[]);

metric = 'bic';
idx = 1:11;     % Humans only

m1 = modelsummary{1}.(metric)(idx,1);
m2 = modelsummary{2}.(metric)(idx,1:2);
m3 = modelsummary{3}.(metric)(idx,1:2);

msep = bsxfun(@plus, m1, m2 + m3);
mjoint = modelsummary{4}.(metric)(idx,1:2);

[msep,mjoint]




end
