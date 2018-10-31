function [bms,mfits] = VestBMS_modelComparison(metric,type,mfits,metricname,priorstrength,flags)
%VESTBMS_MODELCOMPARISON Main model comparison.

if nargin < 1 || isempty(metric); metric = 'aicc'; end
if nargin < 2 || isempty(type); type = 1; end
if nargin < 3 || isempty(mfits); mfits = []; end
if nargin < 4 || isempty(metricname); metricname = metric; end
if nargin < 5 || isempty(priorstrength); priorstrength = 'factor'; end
if nargin < 6 || isempty(flags); flags = [1 1 1 0 0]; end
if numel(flags) < 5; flags(5) = 0; end

baypmflag = flags(5);   % BAY probability matching instead of BAY model average

if isempty(mfits)
    mfits = load('VestBMS_modelfits');
end

bms = [];
metricString = [' (' metricname ')'];

switch type
    case 0
        h = plotify([1 2 3],'margins',[0.05 0.05,0.1 0.05],'labels',{'a','b','c'});
        titleaxes = [1 2 3];
        do_bms = 0;
    case 1
        h = plotify([repmat([1 1 1, 2 2 2, 3 3 3],3,1); [4 5 6, 7 8 9, 10 11 12]],'gutter',[0.05,0.1],'margins',[0.05 0.05,0.1 0.05],'labels',{'a','b','c'});
        titleaxes = [1 2 3];
        do_bms = 1;
    case 2
        nrows = 5; ncols = 3; padding = 5;
        grid = [makepanel(nrows,ncols,padding), zeros(nrows,1), makepanel(nrows,ncols,padding)+15, zeros(nrows,1), makepanel(nrows,ncols,padding)+30];        
        h = plotify(grid,'gutter',[0.03,0.03],'margins',[0.05 0.05,0.1 0.05]);
        titleaxes = [2 17 32];
        do_bms = 2;
end

%% BISENSORY ESTIMATION TASK
if flags(1)
    if baypmflag == 2
        modelnames = {'BPD-C','BPMD-C','FFD-C','CXD-C','BPD','BPMD','FFD','CXD','BP-C','BPM-C','FF-C','CX-C','BP','BPM','FF','CX'};
        M = numel(modelnames); priorweight = [0.5 0.5 1 1, 0.5 0.5 1 1, 0.5 0.5 1 1, 0.5 0.5 1 1];
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 1 1 0 0 0 0 1 1 1 1 0 0 0 0; 0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1];
        factornames{2} = {'Bayes','Fixed','Fusion'};
        factors{2} = [1 1 0 0, 1 1 0 0, 1 1 0 0, 1 1 0 0; 0 0 0 1, 0 0 0 1, 0 0 0 1, 0 0 0 1; 0 0 1 0, 0 0 1 0, 0 0 1 0, 0 0 1 0];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1 1 1 1 1, 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0, 1 1 1 1 1 1 1 1];
    else
        if baypmflag == 1
            modelnames = {'BPMD-C','FFD-C','CXD-C','BPMD','FFD','CXD','BPM-C','FF-C','CX-C','BPM','FF','CX'};
        else
            modelnames = {'BPD-C','FFD-C','CXD-C','BPD','FFD','CXD','BP-C','FF-C','CX-C','BP','FF','CX'};            
        end
        M = numel(modelnames); priorweight = ones(1,M);
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 1 0 0 0 1 1 1 0 0 0; 0 0 0 1 1 1 0 0 0 1 1 1];
        factornames{2} = {'Bayes','Fixed','Fusion'};
        factors{2} = [1 0 0, 1 0 0, 1 0 0, 1 0 0; 0 0 1, 0 0 1, 0 0 1, 0 0 1; 0 1 0, 0 1 0, 0 1 0, 0 1 0];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1 1 1, 0 0 0 0 0 0; 0 0 0 0 0 0 1 1 1 1 1 1];
    end
    switch type
        case 1
            masks = [];
            hg = [h(1),h(4:6)];
            factorfixed = [];
        case 2
            masks = {ones(1,M),[1 1 1 0 0 0 1 1 1 0 0 0],[0 0 0 1 1 1 0 0 0 1 1 1],[1 1 1 1 1 1, 0 0 0 0 0 0],[0 0 0 0 0 0 1 1 1 1 1 1]};
            hg = {[0,h(1:3)], [0,h(4:6)], [0,h(7:9)], [0,h(10:12)], [0,h(13:15)]};
            factorfixed = [0 0 0; 1 0 0; 2 0 0; 0 0 1; 0 0 2];
    end
    
    if do_bms
        [bms{1},fac{1}] = ModelWork_factorialComparison(mfits.modelsummary_biml,metric,modelnames,'BPD',priorweight,factors,factornames,hg,masks,factorfixed,priorstrength);
        bms{1}.factor = fac{1};
        axes(h(titleaxes(1))); title(['Bisensory estimation only', metricString]);
    else
        axes(h(1));
        ModelPlot_compare(mfits.modelsummary_biml,metric,'BPD','mean',[],modelnames);
    end
    drawnow;
    clear factors factornames masks;
end

%% UNITY JUDGMENT TASK
if flags(2)
    if baypmflag == 2
        modelnames = {'BPD-C','BPPD-C','CX-C','BPD','BPPD','CX','BP-C','BPP-C','BP','BPP','FF'};
        M = numel(modelnames); priorweight = [0.5 0.5 2 0.5 0.5 2, 0.5 0.5 0.5 0.5 4];
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 1 0 0 0 1 1 0 0 1; 0 0 0 1 1 1 0 0 1 1 1];
        factornames{2} = {'Bayes','Fixed','Fusion'};
        factors{2} = [1 1 0 1 1 0 1 1 1 1 0; 0 0 1 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1 1 1 0 0 0 0 1; 0 0 1 0 0 1 1 1 1 1 1];
    else
        if baypmflag == 1
            modelnames = {'BPPD-C','CX-C','BPPD','CX','BPP-C','BPP','FF'};
        else
            modelnames = {'BPD-C','CX-C','BPD','CX','BP-C','BP','FF'};            
        end
        M = numel(modelnames); priorweight = [1 2 1 2, 1 1 4];
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 0 0 1 0 1; 0 0 1 1 0 1 1];
        factornames{2} = {'Bayes','Fixed','Fusion'};
        factors{2} = [1 0 1 0 1 1 0; 0 1 0 1 0 0 0; 0 0 0 0 0 0 1];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1 0 0 1; 0 1 0 1 1 1 1];
    end

    switch type
        case 1
            masks = [];
            hg = [h(2),h(7:9)];
            % priorweight = ones(1,7);
            factorfixed = [];
        case 2
            masks = {ones(1,M),[1 1 0 0 1 0 1],[0 0 1 1 0 1 1],[1 1 1 1 0 0 1],[0 1 0 1 1 1 1]};
            hg = {[0,h(16:18)], [0,h(19:21)], [0,h(22:24)], [0,h(25:27)], [0,h(28:30)]};
            priorweight = [1 2 1 2, 1 1 4; 1 2 0 0, 1 0 2; 0 0 1 2, 0 1 2; 1 1 1 1, 0 0 2; 0 1 0 1, 1 1 2];
            factorfixed = [0 0 0; 1 0 0; 2 0 0; 0 0 1; 0 0 2];
    end

    if do_bms
        [bms{2},fac{2}] = ModelWork_factorialComparison(mfits.modelsummary_bimu,metric,modelnames,'BPD',priorweight,factors,factornames,hg,masks,factorfixed,priorstrength);
        bms{2}.factor = fac{2};
        axes(h(titleaxes(2))); title(['Unity judgements only', metricString]);
    else
        axes(h(2));
        ModelPlot_compare(mfits.modelsummary_bimu,metric,'BPD','mean',[],modelnames);
    end
    drawnow;
    clear factors factornames;
end

%% JOINT FIT
if flags(3)
    if baypmflag == 2
        modelnames = {'BPD-C','BPMD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPD','BPMD','CXD','BPFDs','CXFDs','BP-C','BPM-C','BPF-Cs','CX-C','CXF-Cs','BP','BPM','BPFs','CX','CXFs'};        
        M = numel(modelnames); priorweight = [0.5 0.5 1 1 1, 0.5 0.5 1 1 1, 0.5 0.5 1 1 1, 0.5 0.5 1 1 1];
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 1 1 1, 0 0 0 0 0, 1 1 1 1 1, 0 0 0 0 0; 0 0 0 0 0, 1 1 1 1 1, 0 0 0 0 0, 1 1 1 1 1];
        factornames{2} = {'Bayes','Fixed','Bayes/Fusion','Fixed/Fusion'};
        factors{2} = [1 1 0 0 0, 1 1 0 0 0, 1 1 0 0 0, 1 1 0 0 0; 0 0 1 0 0, 0 0 1 0 0, 0 0 1 0 0, 0 0 1 0 0; 0 0 0 1 0, 0 0 0 1 0, 0 0 0 1 0, 0 0 0 1 0; 0 0 0 0 1, 0 0 0 0 1, 0 0 0 0 1, 0 0 0 0 1];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1 1, 1 1 1 1 1, 0 0 0 0 0, 0 0 0 0 0; 0 0 0 0 0, 0 0 0 0 0, 1 1 1 1 1, 1 1 1 1 1];
    else
        if baypmflag == 1
            modelnames = {'BPMD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPMD','CXD','BPFDs','CXFDs','BPM-C','BPF-Cs','CX-C','CXF-Cs','BPM','BPFs','CX','CXFs'};                    
        else
            modelnames = {'BPD-C','CXD-C','BPFD-Cs','CXFD-Cs','BPD','CXD','BPFDs','CXFDs','BP-C','BPF-Cs','CX-C','CXF-Cs','BP','BPFs','CX','CXFs'};
        end
        M = numel(modelnames); priorweight = ones(1,M);
        factornames{1} = {'Constant','Eccentric'};
        factors{1} = [1 1 1 1, 0 0 0 0, 1 1 1 1, 0 0 0 0; 0 0 0 0, 1 1 1 1, 0 0 0 0, 1 1 1 1];
        factornames{2} = {'Bayes','Fixed','Bayes/Fusion','Fixed/Fusion'};
        factors{2} = [1 0 0 0, 1 0 0 0, 1 0 0 0, 1 0 0 0; 0 1 0 0, 0 1 0 0, 0 1 0 0, 0 1 0 0; 0 0 1 0, 0 0 1 0, 0 0 1 0, 0 0 1 0; 0 0 0 1, 0 0 0 1, 0 0 0 1, 0 0 0 1];
        factornames{3} = {'Empirical','Independent'};
        factors{3} = [1 1 1 1, 1 1 1 1, 0 0 0 0, 0 0 0 0; 0 0 0 0, 0 0 0 0, 1 1 1 1, 1 1 1 1];
    end
    switch type
        case 1
            masks = [];
            factorfixed = [];
            hg = [h(3),h(10:12)];
        case 2
            masks = {ones(1,M),[1 1 1 1, 0 0 0 0, 1 1 1 1, 0 0 0 0],[0 0 0 0, 1 1 1 1, 0 0 0 0, 1 1 1 1],[1 1 1 1, 1 1 1 1, 0 0 0 0, 0 0 0 0],[0 0 0 0, 0 0 0 0, 1 1 1 1, 1 1 1 1]};
            factorfixed = [0 0 0; 1 0 0; 2 0 0; 0 0 1; 0 0 2];
            hg = {[0,h(31:33)], [0,h(34:36)], [0,h(37:39)], [0,h(40:42)], [0,h(43:45)]};
    end

    if do_bms
        [bms{3},fac{3}] = ModelWork_factorialComparison(mfits.modelsummary_joint,metric,modelnames,'BPD',priorweight,factors,factornames,hg,masks,factorfixed,priorstrength);
        bms{3}.factor = fac{3};
        axes(h(titleaxes(3))); title(['Joint fits', metricString]);
    else
        axes(h(3));
        ModelPlot_compare(mfits.modelsummary_joint,metric,'BPD','mean',[],modelnames);
    end
    drawnow;
    clear factors factornames;
end

%% UNISENSORY MEASUREMENT TASK
if flags(4)
    modelnames = {'BP-C','BP'}; 
    M = numel(modelnames); priorweight = ones(1,M);
    factornames{1} = {'Constant','Eccentric'};
    factors{1} = [1 0; 0 1];
    switch type
        case 1
            masks = [];
            hg = [];
        case 2
            masks = [];
            hg = [];
    end
    
    if do_bms
        [bms{4},fac{4}] = ModelWork_factorialComparison(mfits.modelsummary_uni,metric,modelnames,'BP',priorweight,factors,factornames,hg,masks,[],priorstrength);
        bms{4}.factor = fac{4};
        axes(h(titleaxes(1))); title(['Unisensory measurement', metricString]);
    else
        axes(h(1));
        ModelPlot_compare(mfits.modelsummary_uni,metric,'BP','mean',[],modelnames);
    end
    drawnow;
    clear factors factornames masks;
end


end

%--------------------------------------------------------------------------
function p = makepanel(nrows,ncols,pad)
    p = [];
    n = 1;
    for i = 1:nrows
        row = [];
        for j = 1:ncols
            row = [row, n*ones(1,pad)];
            n = n + 1;
        end    
        p = [p; row];
    end
end