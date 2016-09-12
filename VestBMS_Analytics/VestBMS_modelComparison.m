%VESTBMS_MODELCOMPARISON Main model comparison.
function VestBMS_modelComparison(metric,BMS)

if nargin < 1 || isempty(metric); metric = 'aicc'; end
if nargin < 2 || isempty(BMS); BMS = 0; end

load('VestBMS_modelfits');

if BMS == 2
    h = plotify([repmat([1 1 1, 2 2 2, 3 3 3],3,1); [4 5 6, 7 8 9, 10 11 12]],'gutter',[0.05,0.1],'margins',[0.05 0.05,0.1 0.05],'labels',{'a','b','c'});
else
    h = plotify([1 2 3],'margins',[0.05 0.05,0.1 0.05],'labels',{'a','b','c'});    
end

%% BISENSORY ESTIMATION TASK

axes(h(1));
priorweight = [];
modelnames = {'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD','BP-C','FF-C','CX-C','BP','FF','CX'}; 
priorweight = 1; M = numel(modelnames);
bms = plotcomparison(modelsummary_biml,metric,modelnames,M,'BPD',priorweight,BMS>0);
title('Bisensory estimation only');

if BMS == 2
    axes(h(4));
    factors = [1 1 1 0 0 0 1 1 1 0 0 0; 0 0 0 1 1 1 0 0 0 1 1 1];
    plotcomparison(modelsummary_biml,metric,modelnames,M,'BPD',priorweight,bms,{'Const','Ecc'},factors);
    
    axes(h(5));
    factors = [1 0 0, 1 0 0, 1 0 0, 1 0 0; 0 0 1, 0 0 1, 0 0 1, 0 0 1; 0 1 0, 0 1 0, 0 1 0, 0 1 0];
    % modelnames = {{'BPD-C','BP-C','BPD','BP'},{'CXD-C','CX-C','CXD','CX'},{'FFD-C','FFD','FF-C','FF'}};
    plotcomparison(modelsummary_biml,metric,modelnames,M,'BPD',priorweight,bms,{'Bayes','Fixed','Fusion'},factors);
    
    axes(h(6));
    factors = [1 1 1 1 1 1, 0 0 0 0 0 0; 0 0 0 0 0 0 1 1 1 1 1 1];
    % modelnames = {{'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD'},{'BP-C','FF-C','CX-C','BP','FF','CX'}};
    plotcomparison(modelsummary_biml,metric,modelnames,M,'BPD',priorweight,bms,{'Corr','Uncorr'},factors);
end

%% UNITY JUDGMENT TASK

axes(h(2));
priorweight = [];
% modelnames = {'BPd-C','BPPd-C','BPD-C','BPPD-C','BPd','BPPd','BPD','BPPD','FF','CX-C','CX'};
% modelnames = {'BPD-C','BPPD-C','CX-C','BPD','BPPD','FF','CX'};
% modelnames = {'BPD-C','CX-C','BPD','FF','CX'};
modelnames = {'BPD-C','CX-C','BPD','CX','BP-C','BP','FF'}; 
priorweight = [1 2 1 2, 1 1 4]; M = numel(modelnames);
bms = plotcomparison(modelsummary_bimu,metric,modelnames,M,'BPD',priorweight,BMS>0);
title('Unity judgements only');

if BMS == 2
    axes(h(7));
    factors = [1 1 0 0 1 0 1; 0 0 1 1 0 1 1];
    plotcomparison(modelsummary_bimu,metric,modelnames,M,'BPD',priorweight,bms,{'Const','Ecc'},factors);
    
    axes(h(8));
    factors = [1 0 1 0 1 1 0; 0 1 0 1 0 0 0; 0 0 0 0 0 0 1];
    % modelnames = {{'BPD-C','BP-C','BPD','BP'},{'CXD-C','CX-C','CXD','CX'},{'FFD-C','FFD','FF-C','FF'}};
    plotcomparison(modelsummary_bimu,metric,modelnames,M,'BPD',priorweight,bms,{'Bayes','Fixed','Fusion'},factors);
    
    axes(h(9));
    factors = [1 1 1 1 0 0 1; 0 1 0 1 1 1 1];
    % modelnames = {{'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD'},{'BP-C','FF-C','CX-C','BP','FF','CX'}};
    plotcomparison(modelsummary_bimu,metric,modelnames,M,'BPD',priorweight,bms,{'Corr','Uncorr'},factors);
end

%% JOINT FIT

axes(h(3));
priorweight = [];
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs','BPMs'};
% modelnames = {'BPD','BPDs','CXD','CXDs','BPFDs','CXFDs','CXs'};
% modelnames = {'BPD','CXD','BPFDs','CXFDs'};
modelnames = {'BPD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPD','CXD','BPFDs','CXFDs','BP-C','BPF-Cs','CX-C','CXF-Cs','BP','BPFs','CX','CXFs'};
M = numel(modelnames);
bms = plotcomparison(modelsummary_joint,metric,modelnames,M,'BPD',priorweight,BMS>0);
title('Joint fits');

if BMS == 2
    axes(h(10));
    factors = [1 1 1 1, 0 0 0 0, 1 1 1 1, 0 0 0 0; 0 0 0 0, 1 1 1 1, 0 0 0 0, 1 1 1 1];
    %modelnames = {{'BPD-C','CXD-C','BPFD-Cs','CXFD-Cs','BP-C','BPF-Cs','CX-C','CXF-Cs'},{'BPD','CXD','BPFDs','CXFDs','BP','BPFs','CX','CXFs'}};
    plotcomparison(modelsummary_joint,metric,modelnames,M,'BPD',priorweight,bms,{'Const','Ecc'},factors);
    
    axes(h(11));
    factors = [1 0 0 0, 1 0 0 0, 1 0 0 0, 1 0 0 0; 0 1 0 0, 0 1 0 0, 0 1 0 0, 0 1 0 0; 0 0 1 0, 0 0 1 0, 0 0 1 0, 0 0 1 0; 0 0 0 1, 0 0 0 1, 0 0 0 1, 0 0 0 1];
    %modelnames = {{'BPD-C','BPD','BP-C','BP'},{'CXD-C','CXD','CX-C','CX'},{'BPFD-Cs','BPFDs','BPF-Cs','BPFs'},{'CXFD-Cs','CXFDs','CXF-Cs','CXFs'}};
    plotcomparison(modelsummary_joint,metric,modelnames,M,'BPD',priorweight,bms,{'Bayes','Fixed','B-Fusion','F-Fusion'},factors);
    
    axes(h(12));
    factors = [1 1 1 1, 1 1 1 1, 0 0 0 0, 0 0 0 0; 0 0 0 0, 0 0 0 0, 1 1 1 1, 1 1 1 1];
    %modelnames = {{'BPD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPD','CXD','BPFDs','CXFDs'},{'BP-C','BPF-Cs','CX-C','CXF-Cs','BP','BPFs','CX','CXFs'}};
    plotcomparison(modelsummary_joint,metric,modelnames,M,'BPD',priorweight,bms,{'Corr','Uncorr'},factors);
end



end

%--------------------------------------------------------------------------
function bms = plotcomparison(modelsummary,metric,modelnames,M,bestmodel,priorweight,BMS,groupnames,factors)

if nargin < 6 || isempty(priorweight); priorweight = 1; end
if nargin < 7 || isempty(BMS); BMS = 0; end
if nargin < 8; groupnames = []; end
if nargin < 9; factors = []; end

if (isnumeric(BMS) || islogical(BMS)) && BMS == 0
    ModelPlot_compare(modelsummary,metric,bestmodel,'mean',[],modelnames);
    bms = [];
else        % Bayesian Model Selection for group studies
    alpha0 = priorweight./sum(priorweight); %/ sqrt(M);
    if (isnumeric(BMS) || islogical(BMS)) && BMS == 1; plottype = 'checker'; else plottype = 'factors'; alpha0 = alpha0; end
    modelsummary.groupnames = groupnames;
    [models,tab,bms] = ModelWork_plotModelComparison(modelsummary,metric,alpha0, ...
        'PlotType',plottype,'ModelList',modelnames,'SSorder',1,'BMSorder',0,'Factors',factors);
    
    if isstruct(BMS)
        ylabel('');
        set(gca,'FontSize',10,'Ytick',[0 0.5 1]);
    else
        bms.modelnames
        bms.exp_r
        h = get(gca,'xlabel');
        set(h,'FontSize',12);
    end
end

end