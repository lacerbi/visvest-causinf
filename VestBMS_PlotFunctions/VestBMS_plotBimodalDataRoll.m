% VESTBMS_PLOTBIMODALDATA plot bimodal data for one or more subjects.
function [fig,gendata] = VestBMS_plotBimodalDataRoll(data,type,mfit,ngen,flags,hg,fontsize,axesfontsize)

if ~exist('mfit', 'var'); mfit = []; end
% Number of generated datasets, per subject
if ~exist('ngen', 'var') || isempty(ngen); ngen = 30; end
if ~exist('flags', 'var') || isempty(flags); flags = 0; end
if ~exist('hg', 'var'); hg = []; end
if ~exist('fontsize', 'var') || isempty(fontsize); fontsize = 14; end
if ~exist('axesfontsize', 'var') || isempty(axesfontsize); axesfontsize = 12; end

plotdata = flags(1);

plots = VestBMS_defaults('plots');  % Get default plot info

% Check dataset
if isstruct(data) && isfield(data, 'X'); data = {data}; end

% Plot average across all provided datasets
subjs = 0;
% subjs = 1:numel(data);

nNoise = 3; % Three levels of noise

fig.prefix = 'VestBMS'; % Program name

if type == 2; nEcc = 5; else nEcc = 3; end

if subjs(1) == 0
    fig.panelgraph = (1:nEcc);
    fig.intborder = [0 0.1]; % Internal border
    nRows = 1;
else
    fig.panelgraph = reshape(1:nEcc*length(subjs),length(subjs),nEcc)';
    fig.intborder = [0.01 0.015]; % Internal border
    nRows = length(subjs);
    plotdata = 0;
end

% stringnoise = {'|$\bar{s}$| = 0°, 5°', '|s| = 10°, 15°', '|s| = 20°, 25°', 'Vestibular'};

yLim = [0 1; 0 1; 0 1];

if type == 2
    % xxx = [-45,-40,-35,-30,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,45];    
    xxx = [-32.5,-27.5,-25,-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,32.5];
    xtick = [-30,-15,0,15,30];
    xticklabel = {'-30°','-15°','0°','15°','30°'};
    xstring = '$s_{vest}$';
    xinterpreter = 'LaTeX';
    % xticklabel = {'-30°','-20°','-10°','0°','10°','20°','30°'};
    ylabels = {'Fraction reported ''right''', 'SD'};
    stringecc = {'$s_{vis} \le -20^\circ$','$-20^\circ < s_{vis} < -5^\circ$','$-5^\circ \le s_{vis} \le 5^\circ$','$5^\circ < s_{vis} < 20^\circ$','$s_{vis} \ge 20^\circ$'};
    xLim = [-35 35; -35 35; -35 35];
    nDisparity = 25;
    legendloc = 'SouthEast';
else
    xxx = [-40,-20,-10,-5,0,5,10,20,40];
    xtick = [-40,-20,-10,-5,0,5,10,20,40];
    xstring = 'Disparity';
    xinterpreter = 'TeX';
    xticklabel = {'-40°','-20°','-10°','-5°','0°','5°','10°','20°','40°'};
    ylabels = {'Fraction reported ''unity''', 'SD'};
    stringecc = {'$\left|\bar{s}\right| = 0^\circ, 5^\circ$', '$\left|\bar{s}\right| = 10^\circ, 15^\circ$', '$\left|\bar{s}\right| = 20^\circ, 25^\circ$', 'Vestibular'};
    xLim = [-45 45; -45 45; -45 45];
    nDisparity = 9;
    legendloc = 'South';
end
binfuns = {'@(y) nanmean(y)', '@(y) nanstd(y)'};

fig.panels = [];

for iEcc = 1:nEcc
    panel = []; iPlot = 1;
    
    for iRow = 1:nRows        
        nid = subjs(iRow);
        if subjs(1) > 0; panel = []; iPlot = 1; end

        for iNoise = 1:nNoise
            % panel = []; iPlot = 1;

            % xsource = ['1:size(thisdata.X.bimflatsymm{' num2str(iNoise) '}{' num2str(type) '},2)'];
            if type == 2
                xsource = ['[' num2str(xxx) ']'];
                ysource = ['thisdata.X.bimright{' num2str(iNoise) '}{' num2str(type) '}((1:' num2str(nDisparity) ') + ' num2str(nDisparity*(iEcc-1)) ')'];
            else
                xsource = ['[' num2str(xxx) ']'];
                ysource = ['thisdata.X.bimflatsymm{' num2str(iNoise) '}{' num2str(type) '}((1:' num2str(nDisparity) ') + ' num2str(nDisparity*(iEcc-1)) ')'];
            end

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
            %if type == 3
            %    panel.plots{iPlot}.source.zfun = '@(x,y,z) 2-z'; 
            %end
            panel.plots{iPlot}.color = plots.NoiseColors(iNoise,:);
            panel.plots{iPlot}.linecolor = plots.NoiseColors(iNoise,:);
            panel.plots{iPlot}.linewidth = 2;
            % panel.plots{iPlot}.errorwidth = 0;
            panel.plots{iPlot}.type = 'errorbar';
            panel.plots{iPlot}.binshift = 0.3*(iNoise-2);

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
            panel.yTick = [0 1];
            if iEcc == 1
                panel.yTickLabel = {'0','1'};            
                % panel.yTickLabel = {'-45°','-30°','-15°','0°','15°','30°','45°'};
            else
                panel.yTickLabel = {'',''};            
            end
            panel.plotzero = 0;

            if iRow == 1; panel.title = stringecc{iEcc}; end
            if iRow == nRows; panel.xlabel = xstring; panel.xlabelinterpreter = xinterpreter; end
            if iEcc == 1 && iRow == nRows; panel.ylabel = ylabels{1}; end

            % Add model fit to panel
            if ~isempty(mfit)
                %for jPlot = 1:length(panel.plots)
                    thisplot = panel.plots{iPlot};
                    iPlot = iPlot + 1;
                    if ~isfield(thisplot.source, 'method')
                        copyplot = thisplot;
                        copyplot.source.type = 'model';
                        % copyplot.interp = 1;
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
        if subjs(1) > 0; fig.panels{end+1} = panel; end % Add panel to figure
        
    end
    if subjs(1) == 0; fig.panels{end+1} = panel; end; % Add panel to figure
end

% fig.hg = hg
[fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen);

for iPanel = 1:numel(fig.hg)
    axes(fig.hg(iPanel));
    title(stringecc{iPanel},'Interpreter','LaTeX');
%    plot([0,0], ylim, 'k:', 'LineWidth', 1);
end

% axes(fig.hg(end));
% for iNoise = 1:nNoise
%    col = 0.5 + 0.5*plots.NoiseColors(iNoise,:);
%    hpatch(iNoise) = patch([0 0], [0 0], col);
%end
%hl = legend(hpatch,'High coherence','Medium coherence','Low coherence');
%set(hl,'Location',legendloc,'box','off');


% h = text(0, 0, 'Probability of rightward response');
% h = text(0, 0, 'Probability of rightward response (monkey subjects)');
%h = text(0, 0, 'Probability of rightward response (human subjects)');
%set(h,'HorizontalAlignment','Center','Position',[0,70],'HandleVisibility','off','FontSize',fontsize);

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