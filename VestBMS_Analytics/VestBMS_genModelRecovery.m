%VESTBMS_GENMODELRECOVERY
function [data,mbag,pcommon,kcommon] = VestBMS_genModelRecovery(task,N,mbag)
%VESTBMS_GENMODELRECOVERY Generate fake datasets for model recovery

% Number of fake datasets per subject
if nargin < 2 || isempty(N); N = 10; end

% Pass external MBAG
if nargin < 3; mbag = []; end

switch lower(task(1))
    case 'u'    % Generate datasets for unity judgement task
        
        % Collect model fits
        if isempty(mbag)
            mbag = ModelWork_collectFits('VestBMS','bisensory-u*',[],[]);
        end
        
        % We consider four models: 
        % ecc Bayesian (BP), ecc fixed (CX), fusion (FF), best constant (BP-C)
        modelnames = {'BP','CX','FF','BP-C'};
        
        filename = 'VestBMS_fakedata_unity';
                                
    case 'l'    % Generate datasets for bisensory localization task
        
        % Collect model fits
        if isempty(mbag)
            mbag = ModelWork_collectFits('VestBMS','bisensory-l*',[],[]);
        end
        
        % We consider four models: 
        % ecc Bayesian probabilistic (BPM), ecc fixed (CX), fusion (FF), best Bayesian constant (BPM-C)
        modelnames = {'BPM','CX','FF','BPM-C'};
        
        filename = 'VestBMS_fakedata_bimloc';

        [pcommon,kcommon] = getcommonfromunity();
        
        phat = betafit(pcommon);
        pcommon_alpha = phat(1);
        pcommon_beta = phat(2);
        
        kcommon_mu = mean(kcommon);
        kcommon_sigma = std(kcommon);
        
        for i = 1:numel(mbag.bag)
            m = mbag.bag{i};
            
            % Randomly generate p_common and replace it
            idx = find(strcmp(m.mp.params,'pcommon'),1);
            if ~isempty(idx)
                m.maptheta(idx) = betarnd(pcommon_alpha,pcommon_beta);
                if ~isempty(m.sampling) && ~isempty(m.sampling.samples)
                    m.sampling.samples(:,idx) = betarnd(pcommon_alpha,pcommon_beta,[size(m.sampling.samples,1),1]);
                end
            end
            
            % Randomly generate k_common and replace it
            idx = find(strcmp(m.mp.params,'kcommon'),1);
            if ~isempty(idx)
                m.maptheta(idx) = normrnd(kcommon_mu,kcommon_sigma);
                if ~isempty(m.sampling) && ~isempty(m.sampling.samples)
                    m.sampling.samples(:,idx) = normrnd(kcommon_mu,kcommon_sigma,[size(m.sampling.samples,1),1]);
                end
            end
            
            mbag.bag{i} = m;
        end
        
    case 'j'    % Generate datasets for joint fits
        
        % We consider six models:
        % (best four model according to BMS with LOO, two variants of the
        %  best model wrt sensory noise and prior)
        modelnames = {'CXD','BPD','BPFs','CXF-Cs','CX','CXD-C'};
        filename = 'VestBMS_fakedata_joint';
        
        
    otherwise
        error('Unknown task. Specify either (U)nity judgment, bimodal (L)ocalization, or (J)oint fits.');
end

% data = [];
data = ModelWork_genModelRecovery(mbag,N,modelnames);

save([filename '.mat'],'data');

end

%--------------------------------------------------------------------------

function [pcommon,kcommon] = getcommonfromunity()
%GETCOMMONFROMUNITY Get p_common and k_common from unity fits

mbag = ModelWork_collectFits('VestBMS','bisensory-u*',[],[]);
modelsummary = ModelWork_summary(mbag);

% modelnames = {'BP','BPP','CX','FF','BPP-C'};

%for i = 1:numel(modelnames)
%    models(i,:) = modelsummary.models(strcmp(modelsummary.modelnames,modelnames{i}),:);
%end
% Use all models
models = modelsummary.models;

% Get all datasets
mfit = ModelBag_get(mbag,modelsummary.dataid(:,:),models,modelsummary.cnd);

% Get p_common and k_common parameters for each model
pcommon = []; kcommon = [];
for i = 1:numel(mfit)
    params = mfit{i}.mp.params;
    idx = find(strcmp(params,'pcommon'),1);
    if ~isempty(idx); pcommon = [pcommon, mfit{i}.maptheta(idx)]; end
    idx = find(strcmp(params,'kcommon'),1);
    if ~isempty(idx); kcommon = [kcommon, mfit{i}.maptheta(idx)]; end
end

end