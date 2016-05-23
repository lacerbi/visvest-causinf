% VESTBMS_PLOTUNIMODALDATANEW plot unimodal data for one or more subjects.
function [fig,gendata] = VestBMS_plotUnimodalData(data,mfit,ngen,flags,fontsize,axesfontsize)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = [1 1 0 0]; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 14; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 12; end

plotdata = flags(1);
plotdots = flags(2);
plotbimodal = flags(3);
% robustfit = flags(5);

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

nid = 0; % Plot average across subjects

if plotbimodal; nNoise = 3; else nNoise = 4; end

subjs = 0;

fig.prefix = 'VestBMS'; % Program name

leftNoise = 1;
bottmNoise = 1:nNoise;

if subjs(1) == 0
    if nNoise == 4
        fig.panelgraph = [1 2; 3 4];
        nRows = 1;
        leftNoise = [1 3];
        bottomNoise = [3 4];
    else
        fig.panelgraph = 1:nNoise;
        nRows = 1;
    end
    fig.intborder = [0.05 0.1]; % Internal border noise
else
    fig.panelgraph = reshape(1:nNoise*length(subjs),nNoise,length(subjs))';
    fig.intborder = [0.05 0.015]; % Internal border noise
    nRows = length(subjs);
    plotdata = 0;
end

stringnoise = {'Visual, low noise', 'Visual, med noise', 'Visual, high noise', 'Vestibular'};

xLim = [-25 25; -25 25];
yLim = [-1.2 1.2; 0 1];

ylabels = {'Pr(respond right)', 'SD'};    
binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];
for iRow = 1:nRows
    nid = subjs(iRow);
    
    for iNoise = 1:nNoise
        panel = []; iPlot = 1;
        
        if ~plotbimodal || 1
            xsource = ['thisdata.X.unimodal{' num2str(iNoise) '}(:, 2)'];
            ysource = ['thisdata.X.unimodal{' num2str(iNoise) '}(:, 3)'];
        else
            modality = plotbimodal;
            xsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(abs(thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,3) - thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,4)) < 1e-6, 3)'];
            ysource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(abs(thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,3) - thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,4)) < 1e-6, 5)'];               
        end
        
        if plotdata
            ysourcedots = ['thisdata.X.unimodal{' num2str(iNoise) '}(:, 3)*1.1'];
            panel.plots{iPlot} = newdataplot(xsource,ysourcedots,nid,[]);
            panel.plots{iPlot}.source.method = 'pool';            
            panel.plots{iPlot}.markertype = '.';
            panel.plots{iPlot}.markersize = 3;
            panel.plots{iPlot}.xjitter = 0.5;
            panel.plots{iPlot}.yjitter = 0.02;
            panel.plots{iPlot}.color = 0.5*[1 1 1];
%            panel.plots{iPlot}.color = col(iNoise,:);
            iPlot = iPlot + 1;
        end
        
        % Mean data plot
        panel.plots{iPlot} = newdataplot(xsource,ysource,nid,binfuns{1});        
        if plotdots; panel.plots{iPlot}.type = 'dots'; else panel.plots{iPlot}.type = 'line'; end
        if plotbimodal; panel.plots{iPlot}.binshift = -0.2; end
        if nRows > 1; panel.plots{iPlot}.markersize = 4; end
                
        % Mean bimodal data plot
        if plotbimodal
            iPlot = iPlot + 1;
            modality = plotbimodal;
            xsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(abs(thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,3) - thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,4)) < 1e-6, 3)'];
            ysource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(abs(thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,3) - thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(modality) '}(:,4)) < 1e-6, 5)'];
            panel.plots{iPlot} = newdataplot(xsource,ysource,nid,binfuns{1});            
            panel.plots{iPlot}.binshift = 0.2;
            panel.plots{iPlot}.color = [0.5 0.5 0.5];
        end
        
        % Panel cosmetics        
        panel.fontsize = fontsize;
        panel.axesfontsize = axesfontsize;
        panel.xLim = xLim(1, :); panel.yLim = yLim(1, :);
        panel.xTick = [-20 -10 0 10 20]; 
        if iRow == nRows
            panel.xTickLabel = {'-20','-10','0°','10','20'};
        else
            panel.xTickLabel = {'','',''};            
        end
        panel.yTick = [-1 0 1];
        if any(iNoise == leftNoise)
            panel.yTickLabel = {'0','0.5','1'};
        else
            panel.yTickLabel = {'','',''};            
        end
        % panel.yTick = [-1.15 -1 0 1 1.15]; panel.yTickLabel = {'L','0','0.5','1','R'};
        
        if iRow == 1; panel.title = stringnoise{iNoise}; end
        if iRow == nRows && any(iNoise == bottomNoise); panel.xlabel = 'Motion direction (deg)'; end
        if any(iNoise == leftNoise) && iRow == nRows; panel.ylabel = ylabels{1}; end
                
        % Add model fit to panel
        if ~isempty(mfit)
            for jPlot = 1:length(panel.plots)
                thisplot = panel.plots{jPlot};
                if ~isfield(thisplot.source, 'method')
                    copyplot = thisplot;
                    copyplot.source.type = 'model';
                    copyplot.type = 'line';
                    copyplot.color = NaN;
                    copyplot.errColor = 0.8*[1 1 1];
                    copyplot.priority = -1;
                    panel.plots{end+1} = copyplot; % Add panel to figure
                end
            end
        end

        fig.panels{end+1} = panel; % Add panel to figure
    end
end

[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen);

return;

%NEWDATAPLOT Basic template for data plot
function newplot = newdataplot(xsource,ysource,nid,binfun)
    newplot.source.type = 'data';
    newplot.source.x = xsource;
    newplot.source.y = ysource;               
    newplot.source.nid = nid;            
    newplot.source.bincenters = [-25:2.5:25];
    % newplot.source.bincenters = [-25:2.5:-2.5,-1.25,-0.625,0.625,1.25,2.5:2.5:25];
    newplot.source.binfun = binfun;
end

end