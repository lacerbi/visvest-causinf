function VestBMS_collectInitialPoints()
%VESTBMS_COLLECTINITIALPOINTS Collect starting points for optimization

erasefields = {'X','mp','infostruct','optimization','sampling','metrics'};

mbag = ModelWork_collectFits('VestBMS','bisensory-unity*',[],[]);
for i = 1:numel(mbag.bag)
    for f = erasefields; mbag.bag{i}.(f{:}) = []; end
end
save('VestBMS_starting_points_unity.mat','mbag');

mbag = ModelWork_collectFits('VestBMS','bisensory-loc*',[],[]);
for i = 1:numel(mbag.bag)
    for f = erasefields; mbag.bag{i}.(f{:}) = []; end
end
save('VestBMS_starting_points_localization.mat','mbag');



end