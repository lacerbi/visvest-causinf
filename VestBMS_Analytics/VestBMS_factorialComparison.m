function mfits = VestBMS_factorialComparison(type,mfits)
%VESTBMS_FACTORIALCOMPARISON Main factorial model comparison.

if nargin < 1 || isempty(type); type = 2; end
if nargin < 2; mfits = []; end

metrics = {'aicc','bic','loocv','marginallike_rlr','marginallike_whmg','marginallike_whmu'};
metricnames = {'AICc','BIC','LOO-CV','MarginalLike_R','MarginalLike_G','MarginalLike_U'};
metricsfile = {'aicc','bic','loocv','mlike-r','mlike-g','mlike-u'};

filetype = {'pdf','png'};

for m = 1:numel(metrics)
    figure;
    [~,mfits] = VestBMS_modelComparison(metrics{m},type,mfits,metricnames{m});
    filename = ['factorial-pxp-' metricsfile{m}];
    
    set(gcf, 'Visible', 'on', 'Position', [1 35 1920 958]);
    for i = 1:numel(filetype)
        export_fig(filetype{i}, [filename '.' filetype{i}], gcf);
    end
end

end
