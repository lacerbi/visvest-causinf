%VESTBMS_GENMODELRECOVERY
function gendata = VestBMS_genModelRecovery(task,N)
%VESTBMS_GENMODELRECOVERY Generate fake datasets for model recovery

% Number of fake datasets per subject
if nargin < 2 || isempty(N); N = 10; end

switch lower(task(1))
    case 'u'    % Generate datasets for unity judgement task
        
        % Collect model fits
        mbag = ModelWork_collectFits('VestBMS','bisensory-u*',[],[]);
        
        % We consider four models: 
        % ecc Bayesian (BP), ecc fixed (CX), fusion (FF), best constant (BP-C)
        modelnames = {'BP','CX','FF','BP-C'};
        
        filename = 'VestBMS_fakedata_unity';
                                
    case 'l'    % Generate datasets for bisensory localization task
        
        % Collect model fits
        mbag = ModelWork_collectFits('VestBMS','bisensory-l*',[],[]);
        
        % We consider four models: 
        % ecc Bayesian (BP), ecc fixed (CX), fusion (FF), best constant (BP-C)
        % modelnames = {'BP','CX','FF','BP-C'};        
        
        filename = 'VestBMS_fakedata_bimloc';
        
    otherwise
        error('Unknown task. Specify either (U)nity judgment or bimodal (L)ocalization.');
end

data = ModelWork_genModelRecovery(mbag,N,modelnames);

save([filename '.mat'],'data');

end