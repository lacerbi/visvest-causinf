%VESTBMS_MODELFITS
function VestBMS_modelFits

[mbag_uni,modelsummary_uni] = ModelWork_collectFits('VestBMS','uni-loc*',[],[]);
[mbag_biml,modelsummary_biml] = ModelWork_collectFits('VestBMS','bim-loc*',[],[]);
[mbag_bimu,modelsummary_bimu] = ModelWork_collectFits('VestBMS','bim-uni*',[],[]);
[mbag_joint,modelsummary_joint] = ModelWork_collectFits('VestBMS','j*',[],[]);

save('VestBMS_modelfits.mat');

end