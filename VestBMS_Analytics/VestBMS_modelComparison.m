%VESTBMS_MODELCOMPARISON Main model comparison.
function VestBMS_modelComparison(metric)

if nargin < 1 || isempty(metric); metric = 'aicc'; end

load('VestBMS_modelfits');

%% BISENSORY ESTIMATION TASK

subplot(1,3,1);
modelnames = {'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD'};
% modelnames = {'BPD-C','BPMD-C','FFD-C','CXD-C','BPTD-C','BPD','BPMD','FFD','CXD'};
%modelnames = {'BPd-C','BPMd-C','BPD-C','BPMD-C','FFd-C','FFD-C','CXd-C','CXD-C','BPTd-C','BPTD-C', ...
%    'BPd','BPMd','BPD','BPMD','FFd','FFD','CXd','CXD'};
ModelPlot_compare(modelsummary_biml,metric,'BPD','mean',[],modelnames);

%% UNITY JUDGMENT TASK

subplot(1,3,2);
% modelnames = {'BPd-C','BPPd-C','BPD-C','BPPD-C','BPd','BPPd','BPD','BPPD','FF','CX-C','CX'};
% modelnames = {'BPD-C','BPPD-C','CX-C','BPD','BPPD','FF','CX'};
modelnames = {'BPD-C','CX-C','BPD','FF','CX'};
ModelPlot_compare(modelsummary_bimu,metric,'BPD','mean',[],modelnames);

%% JOINT FIT

subplot(1,3,3);
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs','BPMs'};
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs'};
modelnames = {'BPD','CXD','BPFDs','CXFDs'};
% modelnames = {'BPDs','CXDs','BPFDs','CXFDs'};
ModelPlot_compare(modelsummary_joint,metric,'BPD','mean',[],modelnames);

end