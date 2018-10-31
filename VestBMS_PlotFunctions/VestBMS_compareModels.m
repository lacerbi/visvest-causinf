%VESTBMS_COMPAREMODELS

if ~exist('modelsummary','var')
    load('modelfits_bim.mat','modelsummary');
end

close all;

action = 'mean';    % Plot average over subjects
baseline = 3;       % Baseline model (FF: forced fusion)

% Human data
subplot(2,2,1);
ModelPlot_compare(modelsummary,'aicc',baseline,action,1:11);
text(-0.3,0.5,'Human','Units','Normalized','FontSize',14);
subplot(2,2,2);
ModelPlot_compare(modelsummary,'bic',baseline,action,1:11);



% Monkey data
subplot(2,2,3);
ModelPlot_compare(modelsummary,'aicc',baseline,action,12:14);
text(-0.3,0.5,'Monkey','Units','Normalized','FontSize',14);
subplot(2,2,4);
ModelPlot_compare(modelsummary,'bic',baseline,action,12:14);