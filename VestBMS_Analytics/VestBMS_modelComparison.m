%VESTBMS_MODELCOMPARISON Main model comparison.
function VestBMS_modelComparison(metric)

if nargin < 1 || isempty(metric); metric = 'aicc'; end

load('VestBMS_modelfits');

%% BISENSORY ESTIMATION TASK

subplot(1,3,1);
priorweight = [];
% modelnames = {'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD'};
modelnames = {'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD','BP-C','FF-C','CX-C','BP','FF','CX'}; priorweight = 1;
% modelnames = {'BPD-C','BPMD-C','FFD-C','CXD-C','BPTD-C','BPD','BPMD','FFD','CXD'};
%modelnames = {'BPd-C','BPMd-C','BPD-C','BPMD-C','FFd-C','FFD-C','CXd-C','CXD-C','BPTd-C','BPTD-C', ...
%    'BPd','BPMd','BPD','BPMD','FFd','FFD','CXd','CXD'};
plotcomparison(modelsummary_biml,metric,modelnames,'BPD',priorweight);

%% UNITY JUDGMENT TASK

subplot(1,3,2); 
priorweight = [];
% modelnames = {'BPd-C','BPPd-C','BPD-C','BPPD-C','BPd','BPPd','BPD','BPPD','FF','CX-C','CX'};
% modelnames = {'BPD-C','BPPD-C','CX-C','BPD','BPPD','FF','CX'};
% modelnames = {'BPD-C','CX-C','BPD','FF','CX'};
modelnames = {'BPD-C','CX-C','BPD','CX','BP-C','BP','FF'}; priorweight = [1 2 1 2, 1 1 4];
plotcomparison(modelsummary_bimu,metric,modelnames,'BPD',priorweight);

%% JOINT FIT

subplot(1,3,3);
priorweight = [];
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs','BPMs'};
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs'};
% modelnames = {'BPD','CXD','BPFDs','CXFDs'};
modelnames = {'BPD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPD','CXD','BPFDs','CXFDs'};
plotcomparison(modelsummary_joint,metric,modelnames,'BPD',priorweight);

end

%--------------------------------------------------------------------------
function plotcomparison(modelsummary,metric,modelnames,bestmodel,priorweight)

if nargin < 5 || isempty(priorweight); priorweight = 1; end

if 0
    ModelPlot_compare(modelsummary,metric,bestmodel,'mean',[],modelnames);
else
    alpha0 = priorweight./sum(priorweight) / numel(modelnames);
    [models,tab,bms] = ModelWork_plotModelComparison(modelsummary,metric,[],modelnames,[],alpha0,1,1)
end

end