% VESTBMS_PLOTBIMODALDISPARITYDATA plot bimodal data for one or more subjects by disparity.
function [fig,gendata] = VestBMS_plotBimodalDisparityData(data,type,mfit,ngen,flags,fontsize,axesfontsize,hg)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = [0 0]; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 18; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 14; end
if ~exist('hg', 'var'); hg = []; end

plotdata = flags(1);
% plotsbar = flags(2);
plotsvis = flags(2);

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

% Plot average across all provided datasets
if numel(data) == 1
    subjs = 1;
else
    subjs = 1:numel(data);
    subjs = 0;
end

nNoise = 3; % Three levels of noise
plots = VestBMS_defaults('plots');  % Get default plot info
col = plots.NoiseColors;            % Colors for different noise levels

% Use provided figure handle
if ~isempty(hg); fig.hg = hg; end

fig.prefix = 'VestBMS'; % Program name

if subjs(1) == 0
    fig.intborder = [0.05 0.1]; % Internal border noise
    nRows = length(type);
    subjs = zeros(1,nRows);
else
    fig.intborder = [0.05 0.015]; % Internal border noise
    nRows = length(subjs);
    type = type*ones(1,nRows);
    plotdata = 0;
end
fig.panelgraph = reshape(1:nNoise*nRows,nNoise,nRows)';

stringnoise = {'Visual, low noise', 'Visual, med noise', 'Visual, high noise', 'Vestibular'};

xLim = [-50 50];
yLim = [0 1; 0 1; 0 1];

ylabels = {'P(vis = right)','P(vest = right)','P(unity = 1)'};    
binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];
for iRow = 1:nRows
    nid = subjs(iRow);
    tt = type(iRow);
    
    for iNoise = 1:nNoise
        panel = []; iPlot = 1;
        
        if plotsvis
            xsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,4)'];            
            % xsource = ['0.5*(thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,3) + thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,4))'];
        else
            xsource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,3) - thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,4)'];
        end
        ysource = ['thisdata.X.bimodal{' num2str(iNoise) '}{' num2str(tt) '}(:,5)'];
                
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
        panel.plots{iPlot} = newdataplot(xsource,ysource,nid,binfuns{1});        
        if plotsvis
            panel.plots{iPlot}.source.bincenters = [-35,-20,-10,-5,0,5,10,20,35];
        end            
        
        if tt == 3
            panel.plots{iPlot}.source.yfun = '@(x,y) 2-y'; 
        else
            panel.plots{iPlot}.source.yfun = '@(x,y) 0.5*(y+1)';
        end
        panel.plots{iPlot}.color = col(iNoise,:);
        panel.plots{iPlot}.type = 'errorbar';
        panel.plots{iPlot}.binshift = 0.25*(iNoise-2);
        panel.plots{iPlot}.linewidth = 4;
                                
        % Panel cosmetics        
        panel.fontsize = fontsize;
        panel.axesfontsize = axesfontsize;
        panel.xLim = xLim(1, :); panel.yLim = yLim(tt, :);
        panel.xTick = [-40 -20 0 20 40];
        % panel.axissquare = 1;
        if iRow == nRows
            panel.xTickLabel = {'-40°','-20°','0°','20°','40°'};
        else
            panel.xTickLabel = {'','','','',''};            
        end
        panel.yTick = [-1 0 1];
        if iNoise == 1
            panel.yTickLabel = {'-1','0','1'};
        else
            panel.yTickLabel = {'','',''};            
        end
        panel.plotzero = 0;
        
        if iRow == 1; panel.title = stringnoise{iNoise}; end
        if iRow == nRows; panel.xlabel = 'Stimulus disparity'; end
        if iNoise == 1; panel.ylabel = ylabels{tt}; end
                
        % Add model fit to panel
        if ~isempty(mfit)
            for jPlot = 1:length(panel.plots)
                thisplot = panel.plots{jPlot};
                if ~isfield(thisplot.source, 'method')
                    copyplot = thisplot;
                    copyplot.source.type = 'model';
                    copyplot.type = 'line';
                    copyplot.color = NaN;
                    % copyplot.errColor = 0.8*[1 1 1];
                    copyplot.errColor = col(iNoise,:);
                    copyplot.priority = -1;
                    copyplot.facealpha = 0.4;
                    copyplot.source.dataids = nid;
                    panel.plots{end+1} = copyplot; % Add panel to figure                    
                end
            end
        end

        fig.panels{end+1} = panel; % Add panel to figure
    end
end

options.flatten = 0;    % Do not compute flattened bins
[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen,options);

plot(xlim,[0.5 0.5],'k--','LineWidth',1);
set(gca,'YTick',0:0.25:1,'YTickLabel',0:0.25:1);


return;

%NEWDATAPLOT Basic template for data plot
function newplot = newdataplot(xsource,ysource,nid,binfun)
    newplot.source.type = 'data';
    newplot.source.x = xsource;
    newplot.source.y = ysource;              
    newplot.source.nid = nid;
    newplot.source.bincenters = [-40:2.5:40];
    % newplot.source.bincenters = [-25:2.5:-2.5,-1.25,-0.625,0.625,1.25,2.5:2.5:25];
    newplot.source.binfun = binfun;
end

end