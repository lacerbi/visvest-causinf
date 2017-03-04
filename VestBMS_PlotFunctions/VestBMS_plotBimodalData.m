% VESTBMS_PLOTBIMODALDATA plot bimodal data for one or more subjects.
function [fig,gendata] = VestBMS_plotBimodalData(data,type,mfit,ngen,flags,hg,fontsize,axesfontsize)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = [0 1]; end
if ~exist('hg', 'var'); hg = []; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 14; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 12; end

plotdata = flags(1);
plot1d = flags(2);

plots = VestBMS_defaults('plots');  % Get default plot info

% Use provided figure handle
if ~isempty(hg); fig.hg = hg; end

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

% Plot average across all provided datasets
subjs = 0;

nNoise = 3; % Three levels of noise

fig.prefix = 'VestBMS'; % Program name

if subjs(1) == 0
    if plot1d
        fig.panelgraph = (1:nNoise)';
    else
        fig.panelgraph = 1:nNoise;
    end
    fig.intborder = [0.05 0.1]; % Internal border noise
    nRows = 1;
else
    fig.panelgraph = reshape(1:nNoise*length(subjs),nNoise,length(subjs))';
    fig.intborder = [0.05 0.015]; % Internal border noise
    nRows = length(subjs);
    plotdata = 0;
end

stringnoise = {'Visual, low noise', 'Visual, med noise', 'Visual, high noise', 'Vestibular'};

if plot1d
    xLim = [0 100; 0 100];
    yLim = [0 1; 0 1];
    zLim = [-1 1; -1 1; 0 1];    
else
    xLim = [-46.4 46.4; -25 25];
    yLim = [-46.4 46.4; -25 25];
    zLim = [-1 1; -1 1; 0 1];
end

ylabels = {'Visual stimulus', 'SD'};    
binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];
for iRow = 1:nRows
    nid = subjs(iRow);
    
    for iNoise = 1:nNoise
        panel = []; iPlot = 1;
        
        xsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(type) '}(:,3)'];
        ysource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(type) '}(:,4)'];
        zsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(type) '}(:,5)'];            
        
        
        if plotdata && 0
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
        panel.plots{iPlot} = newdataplot(xsource,ysource,zsource,nid,binfuns{1});        
        if type == 3
            panel.plots{iPlot}.source.zfun = '@(x,y,z) 2-z'; 
        end
        if plot1d
            panel.plots{iPlot}.type = 'dots2dflat'; 
            panel.plots{iPlot}.color = plots.NoiseColors(iNoise,:);
        else
            panel.plots{iPlot}.color = [1 0 0];
        end
                        
        % Panel cosmetics        
        panel.fontsize = fontsize;
        panel.axesfontsize = axesfontsize;
        panel.plotzero = 0;        
        panel.xLim = xLim(1, :); panel.yLim = yLim(1, :); panel.zLim = zLim(type, :);
        if ~plot1d;  end
        
        if plot1d
            panel.xTick = [1,10:10:90,99];
            panel.yTick = [0 0.5 1];            
            if iRow == nRows; panel.xlabel = 'Stimuli pairs'; end
            panel.ylabel = 'Fraction response ''unity''';            

        else
            panel.axissquare = 1;        
            panel.xTick = [-45 -30 -15 0 15 30 45];
            if iRow == nRows
                panel.xTickLabel = {'-45°','-30°','-15°','0°','15°','30°','45°'};
            else
                panel.xTickLabel = {'','',''};            
            end
        
            panel.yTick = [-45 -30 -15 0 15 30 45];
            if iNoise == 1
                panel.yTickLabel = {'-45°','-30°','-15°','0°','15°','30°','45°'};
            else
                panel.yTickLabel = {'','',''};            
            end
            
            if iRow == 1; panel.title = stringnoise{iNoise}; end
            if iRow == nRows; panel.xlabel = 'Vestibular stimulus'; end
            if iNoise == 1 && iRow == nRows; panel.ylabel = ylabels{1}; end            
        end
                
        % Add model fit to panel
        if ~isempty(mfit)
            for jPlot = 1:length(panel.plots)
                thisplot = panel.plots{jPlot};
                if ~isfield(thisplot.source, 'method')
                    copyplot = thisplot;
                    copyplot.source.type = 'model2d';
                    if plot1d
                        copyplot.type = 'line2dflat'; 
                        copyplot.color = 0.4*plots.NoiseColors(iNoise,:) + 0.6*[1 1 1];  
                        copyplot.linewidth = 4;
                    else
                        copyplot.color = [0 0 0];
                    end
                    % copyplot.color = 'none';
                    panel.plots{end+1} = copyplot; % Add panel to figure
                end
            end
        end

        fig.panels{end+1} = panel; % Add panel to figure
    end
end

[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen);

for iPanel = 1:numel(fig.hg)-1
    axes(fig.hg(iPanel)); hold on;
    if plot1d
        hp = plot(xlim, [0.5 0.5], 'k--', 'LineWidth', 1);
    else
        hp = plot([0,0], ylim, 'k--', 'LineWidth', 1);
    end
    uistack(hp,'bottom');
end

axes(fig.hg(2));

% h = text(0, 0, 'Probability of rightward response');
% h = text(0, 0, 'Probability of rightward response (monkey subjects)');
h = text(0, 0, 'Probability of rightward response (human subjects)');
set(h,'HorizontalAlignment','Center','Position',[0,70],'HandleVisibility','off','FontSize',fontsize);

return;

%NEWDATAPLOT Basic template for data plot
function newplot = newdataplot(xsource,ysource,zsource,nid,binfun)
    newplot.source.type = 'data2d';
    newplot.source.x = xsource;
    newplot.source.y = ysource;              
    newplot.source.z = zsource;              
    newplot.source.nid = nid;
    newplot.source.xbincenters = [-45:2.5:45];
    newplot.source.ybincenters = [-45:2.5:45];
    % newplot.source.bincenters = [-25:2.5:-2.5,-1.25,-0.625,0.625,1.25,2.5:2.5:25];
    newplot.source.binfun = binfun;
end

end