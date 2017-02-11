function mfits = VestBMS_fixModels(mfits)
%VESTBMS_FIXMODELS Add model information.

% Add log flag to log transformed parameters
logparams = {'sigma_vis_low','sigma_vis_med','sigma_vis_high','sigma_vest','kcommon','priorsigma','priorsigmadelta'};

ff = fieldnames(mfits)';
for f = ff
    if ~isfield(mfits.(f{:}),'bag'); continue; end
    % mfits.(f{:}) = ModelWork_loadFields(mfits.(f{:}));
    
    bag = mfits.(f{:}).bag;
    for i = 1:numel(bag)
        if isempty(bag{i}.mp); continue; end
        logflag = zeros(1,numel(bag{i}.mp.params));
        for j = 1:numel(logparams)
            idx = find(strcmp(logparams{j},bag{i}.mp.params),1);
            if isempty(idx); continue; end
            logflag(idx) = 1;
        end
        bag{i}.mp.bounds.logflag = logflag;
    end
    mfits.(f{:}).bag = bag;
end