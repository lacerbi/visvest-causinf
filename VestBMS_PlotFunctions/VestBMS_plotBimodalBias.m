% VESTBMS_PLOTBIMODALPSE plot bimodal data for one or more subjects.
function [fig,gendata] = VestBMS_plotBimodalBias(data,type,mfit,ngen,flags,fontsize,axesfontsize)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = 0; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 14; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 12; end

plotdata = flags(1);

plots = VestBMS_defaults('plots');  % Get default plot info

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

% Plot average across all provided datasets
subjs = 0;

nNoise = 3; % Three levels of noise

fig.prefix = 'VestBMS'; % Program name

if subjs(1) == 0
    fig.panelgraph = 1;
    fig.intborder = [0 0.1]; % Internal border noise
    nRows = 1;
else
    fig.panelgraph = 1:length(subjs);
    fig.intborder = [0.05 0.015]; % Internal border noise
    nRows = length(subjs);
end

yLim = [-20 20];

xxx = [-40,-20,-10:5:10,20,40];

xtick = xxx;
for i = 1:numel(xxx); xticklabel{i} = [num2str(xtick(i)) '°']; end
xstring = 'Disparity';
xinterpreter = '';
ystring = 'Bias';
xLim = [-45 45; -45 45; -45 45];
legendloc = 'NorthWest';

binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];

panel = []; iPlot = 1;

for iRow = 1:nRows
    nid = subjs(iRow);

    for iNoise = 1:nNoise
        % panel = []; iPlot = 1;

        xsource = ['[' num2str(xxx) ']'];
        ysource = ['thisdata.psyrightdelta_mu(' num2str(iNoise) ',:)'];

        % Mean data plot
        panel.plots{iPlot} = newdataplot(xsource,ysource,nid,binfuns{1});        
        %if type == 3
        %    panel.plots{iPlot}.source.zfun = '@(x,y,z) 2-z'; 
        %end
        panel.plots{iPlot}.color = plots.NoiseColors(iNoise,:);
        panel.plots{iPlot}.linecolor = plots.NoiseColors(iNoise,:);
        panel.plots{iPlot}.linewidth = 2;
        panel.plots{iPlot}.errorwidth = 0;
        panel.plots{iPlot}.type = 'errorbar';
        panel.plots{iPlot}.binshift = 0.25*(iNoise-2);

        % Panel cosmetics        
        panel.fontsize = fontsize;
        panel.axesfontsize = axesfontsize;
        panel.xLim = xLim(1, :); panel.yLim = yLim(1, :);
        panel.xTick = xtick;
        % panel.axissquare = 1;
        if iRow == nRows
            panel.xTickLabel = xticklabel;
        else
            panel.xTickLabel = {'','',''};            
        end
        panel.yTick = [-45:15:45];
        panel.yTickLabel = {'-45°','-30°','-15°','0°','15°','30°','45°'};
        panel.plotzero = 1;

        % panel.title = 'PSE';
        if iRow == nRows; panel.xlabel = xstring; panel.xlabelinterpreter = xinterpreter; end
        panel.ylabel = ystring;

        % Add model fit to panel
        if ~isempty(mfit)
            %for jPlot = 1:length(panel.plots)
                thisplot = panel.plots{iPlot};
                iPlot = iPlot + 1;
                if ~isfield(thisplot.source, 'method')
                    copyplot = thisplot;
                    copyplot.source.type = 'model';
                    copyplot.interp = 0;
                    % copyplot.type = 'line';
                    copyplot.type = 'line';
                    copyplot.color = NaN;
                    copyplot.errColor = plots.NoiseColors(iNoise,:);
                    copyplot.priority = -1;
                    copyplot.facealpha = 0.3;
%                    copyplot.color = [0 0 0];
                    % copyplot.color = 'none';
                    panel.plots{end+1} = copyplot; % Add panel to figure
                end
            % end
        end

        iPlot = iPlot + 1;

        % fig.panels{end+1} = panel; % Add panel to figure
    end
end
fig.panels{end+1} = panel; % Add panel to figure

options.psycholeftrightdelta = 1;    % Compute psychometric functions
[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen,options);

title('Vestibular Bias');

% axes(fig.hg(end));
for iNoise = 1:nNoise
    col = 0.5 + 0.5*plots.NoiseColors(iNoise,:);
    hpatch(iNoise) = patch([0 0], [0 0], col);
end
hl = legend(hpatch,'High coherence','Medium coherence','Low coherence');
set(hl,'Location',legendloc,'box','off');

return;

%NEWDATAPLOT Basic template for data plot
function newplot = newdataplot(xsource,ysource,nid,binfun)
    newplot.source.type = 'data';
    newplot.source.x = xsource;
    newplot.source.y = ysource;
    % newplot.source.method = 'pool';

    newplot.source.nid = nid;
    newplot.source.xbincenters = xxx;
    % newplot.source.ybincenters = [1:99];
    newplot.source.bincenters = xxx;
    newplot.source.binfun = binfun;
end

end