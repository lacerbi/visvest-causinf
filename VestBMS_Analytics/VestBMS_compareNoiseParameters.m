function [ratio,ratiosd] = VestBMS_compareNoiseParameters(mbag_uni,mbag_bim)
%VESTBMS_COMPARENOISEPARAMETERS Compare posteriors b/w us and bs sessions.

Nsamples = 1e6;

% Get uni-sensory fits
for i = 1:11; dataid_uni{i} = [i, 8]; end
model_uni = [5 3 1 1 1 1 1 1 1 1 2 1 1 1 1 1];
cnd_uni = [1 2 3 4];
mfit_uni = ModelBag_get(mbag_uni, dataid_uni, model_uni, cnd_uni);
mfit_uni = ModelWork_loadFields('VestBMS',mfit_uni);

% Get bi-sensory fits
dataid_bim = dataid_uni;
model_bim = [5 3 1 1 1 1 1 1 1 1 2 1 1 1 5 1];
cnd_bim = [5 6 7];
mfit_bim = ModelBag_get(mbag_bim, dataid_bim, model_bim, cnd_bim);
mfit_bim = ModelWork_loadFields('VestBMS',mfit_bim);

Nsubjs = numel(mfit_uni);
if numel(mfit_bim) ~= Nsubjs
    error('Mismatch between uni-sensory and bi-sensory model fits.');
end

% params = {'sigma_vis_low','sigma_vis_med','sigma_vis_high','sigma_vest','w_vis','w_vest'};
params = {'sigma_vis_low','sigma_vis_med','sigma_vis_high','sigma_vest'};
Nparams = numel(params);
idx_uni = zeros(1,Nparams);
idx_bim = zeros(1,Nparams);
for i = 1:Nparams
    idx_uni(i) = find(strcmp(params{i},mfit_uni{1}.mp.params),1);
    idx_bim(i) = find(strcmp(params{i},mfit_bim{1}.mp.params),1);
end

ratio = NaN(Nsubjs,Nparams);
ratiosd = NaN(Nsubjs,Nparams);

for i = 1:Nsubjs
    for j = 1:Nparams                
        % Take random set of unisensory and bisensory samples
        su = mfit_uni{i}.sampling.samples(:, idx_uni(j));
        sb = mfit_bim{i}.sampling.samples(:, idx_bim(j));
        su = su(randi(size(su,1),[Nsamples,1]));
        sb = sb(randi(size(sb,1),[Nsamples,1]));        
        if params{j}(1) == 's'
            d = exp(sb-su);
        else
            d = sb./su;
        end
        ratio(i,j) = mean(d);
        ratiosd(i,j) = std(d);        
    end
end

end