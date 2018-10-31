% VESTBMS_PLOTBIMODALPSE plot bimodal data for one or more subjects.
function [fig,gendata] = VestBMS_plotBimodalPSE(data,type,mfit,ngen,flags,fontsize,axesfontsize,hg)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = [0 0]; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 18; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 14; end
if ~exist('hg', 'var'); hg = []; end

plotbias = flags(1);
if numel(flags) < 2; flags(2) = 0; end
plotsd = flags(2);

plots = VestBMS_defaults('plots');  % Get default plot info

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

% Plot average across all provided datasets
subjs = 0;

nNoise = 3; % Three levels of noise

% Use provided figure handle
if ~isempty(hg); fig.hg = hg; end

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

if plotsd
    yLim = [0 30];    
else
    yLim = [-20 20];
end

allStimuli =  {-45:5:-30,-27.5:2.5:-22.5,[-20,-17.5,-15],[-12.5,-10,-7.5],[-5,-2.5],0,[2.5,5],[7.5,10,12.5],[15,17.5,20],22.5:2.5:27.5,30:5:45};
xxx = cellfun(@mean,allStimuli);

% xtick = xxx;
% xticklabel = {'-37.5°','-25°','-17.5°','-10°','-3.75°','0°','3.75°','10°','17.5°','25°','37.5°'};
xtick = [-40,-20,0,20,40];
xticklabel = {'-40°','-20°','0°','20°','40°'};

% xstring = '$\left\langle s_\mathrm{vis} \right\rangle$';
xstring = '$s_\mathrm{vis}$';
xinterpreter = 'LaTeX';
ystring = 'PSE';
%xLim = [-40 40; -40 40; -40 40];
xLim = [-50 50; -50 50; -50 50];
legendloc = 'NorthEast';

binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];

panel = []; iPlot = 1;

for iRow = 1:nRows
    nid = subjs(iRow);

    for iNoise = 1:nNoise
        % panel = []; iPlot = 1;

        xsource = ['[' num2str(xxx) ']'];
        if plotsd
            ysource = ['thisdata.psyright_sigma(' num2str(iNoise) ',:)'];            
        elseif plotbias
            ysource = ['-thisdata.psyright_mu(' num2str(iNoise) ',:)'];            
        else
            ysource = ['thisdata.psyright_mu(' num2str(iNoise) ',:)'];
        end

        % Mean data plot
        panel.plots{iPlot} = newdataplot(xsource,ysource,nid,binfuns{1});        
        %if type == 3
        %    panel.plots{iPlot}.source.zfun = '@(x,y,z) 2-z'; 
        %end
        panel.plots{iPlot}.color = plots.NoiseColors(iNoise,:);
        panel.plots{iPlot}.linecolor = plots.NoiseColors(iNoise,:);
        panel.plots{iPlot}.linewidth = 4;
        panel.plots{iPlot}.errorwidth = [];
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
        panel.plotzero = 2;

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
                    copyplot.facealpha = 0.4;
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


options.psycholeftright = 1;    % Compute psychometric functions
[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen,options);

% title('Vestibular PSE');

% axes(fig.hg(end));
%for iNoise = 1:nNoise
%    col = 0.5 + 0.5*plots.NoiseColors(iNoise,:);
%    hpatch(iNoise) = patch([0 0], [0 0], col);
%end
%hl = legend(hpatch,'High coherence','Medium coherence','Low coherence');
%set(hl,'Location',legendloc,'box','off');

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