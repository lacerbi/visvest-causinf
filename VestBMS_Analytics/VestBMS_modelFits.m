%VESTBMS_MODELFITS
function VestBMS_modelFits

%[mbag_uni,modelsummary_uni] = ModelWork_collectFits('VestBMS','uni-loc*',[],[]);
[mbag_uni,modelsummary_uni] = ModelWork_collectFits('VestBMS','unis*',[],[]);
[mbag_biml,modelsummary_biml] = ModelWork_collectFits('VestBMS','bisensory-l*',[],[]);
[mbag_bimu,modelsummary_bimu] = ModelWork_collectFits('VestBMS','bisensory-u*',[],[]);
[mbag_joint,modelsummary_joint] = ModelWork_collectFits('VestBMS','joint*',[],[]);
[mbag_semi,modelsummary_semi] = ModelWork_collectFits('VestBMS','semi*',[],[]);

save('VestBMS_modelfits.mat');

end