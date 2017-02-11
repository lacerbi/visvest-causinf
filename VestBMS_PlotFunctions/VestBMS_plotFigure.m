function varargout = VestBMS_plotFigure(fig,subfig,mfits,bms,filetype,filesuffix)
%VESTBMS_PLOTFIGURE Plot a specific Figure for the paper or report.
% (Supplementary figures start from 101.)
%
% Usage example:
% CueBMS_plotFigure(fig, data)
% 

fontsize = 18;
axesfontsize = 14;

Nsubjs = 11;

if nargin < 2; subfig = []; end

% Load model fit results
if nargin < 3 || isempty(mfits); mfits = load('VestBMS_modelfits'); end
if nargin < 4; bms = []; end

mfits = VestBMS_fixModels(mfits);
mfits.mbag_uni.legend = 'Unisensory measurement';
mfits.mbag_biml.legend = 'Bisensory implicit CI';
mfits.mbag_bimu.legend = 'Bisensory explicit CI';
mfits.mbag_joint.legend = 'Joint datasets';

% Load datasets
load('VestBMS_data.mat');

plots = VestBMS_defaults('plots');

switch fig
    case {2,3}  % Explicit, implicit inference plot
        
        switch fig
            case 2
                task = 3; stat = 1; bms_idx = 2;
                ystring = 'Fraction responses ''unity''';
                modelnames = {'BPD','BP-C','CX-C','FF'};
                modelnames_text = {'BAY-X-E','BAY-C-I','FIX-C','SFU'};
                mbag = mfits.mbag_bimu;
                legloc = 'NorthWest';   mtxtpos = 0.1;
                factornames = {'Noise','CI Strategy','Prior'};
                flags = [0 0];
            case 3
                task = 2; stat = 2; bms_idx = 1;
                % ystring = 'Fraction responses ''right''';
                modelnames = {'FFD','BPD','CXD-C','BPD-C'};
                modelnames_text = {'FFU-X-E','BAY-X-E','FIX-C-E','BAY-C-E'};
                mbag = mfits.mbag_biml;
                mtxtpos = 0.4;
                factornames = {'Noise','CI Strategy','Prior'};
                flags = 1;  % Plot bias?
                if flags
                    ystring = 'Vestibular bias';
                    legloc = 'SouthEast';
                else
                    ystring = 'Vestibular PSE';
                    legloc = 'SouthWest';
                end
        end
        Ngen = 20;   % Number of generated datasets per subject
        
        % Full datasets for 11 subjects
        if isempty(subfig); subjs = 1:11; else subjs = subfig; end
        data = data(subjs);
                
        % Plot data
        grid = [1 1 1 1 1 1, 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6;  1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 7 7 7 8 8 8; 9 9 9 9 9 9 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8];
        hg = plotify(grid,'Margins',[0.1 0.1],'Labels',{'a','b','','','c'},'FontSize',fontsize);
        hg_a = [hg(1),hg(1),hg(1)];
        if fig == 2
            VestBMS_plotBimodalDisparityData(data,task,[],[],flags,[],[],hg_a);
        elseif fig == 3
            VestBMS_plotBimodalPSE(data,task,[],[],flags,[],[],hg(1));
        end
        title(''); ylabel(ystring);
        for iCol = 1:3; hcol(iCol) = plot(0,0,'-','Color',plots.NoiseColors(iCol,:),'LineWidth',4); end
        hl = legend(hcol,'High reliability','Medium reliability','Low reliability');
        set(hl,'Box','off','FontSize',axesfontsize,'Location',legloc);
        drawnow;
        
        % Plot inset (Figure 3)
        if fig == 3
            bin = 10;
            D.psyright_mu = data{3}.psyright_mu;
            D.psyright_sigma = data{3}.psyright_sigma;
            pos = get(gca,'Position');
            if flags(1)
                D.psyright_mu = -D.psyright_mu;
                newpos = [pos(1) + pos(3)*0.1, pos(2) + pos(4)*0.7, pos(3)*0.325, pos(4)*0.3];
                xb = -25;
            else
                newpos = [pos(1) + pos(3)*0.65, pos(2) + pos(4)*0.7, pos(3)*0.3, pos(4)*0.3];
                xb = [];
            end
            axes('Position',newpos);
            VestBMS_plotPSEdemo(D,bin,xb);
            if flags(1)
                title('Bias = -PSE','FontSize',14);
            end
            % pause
            drawnow;
        end
        
        
        % Plot model fits
        [g,g_modelnames] = VestBMS_stats([1,stat],mfits,bms);
        
        hg_c = hg(5:5+numel(modelnames));
        modelsummary = mfits.modelsummary_bimu;
        for m = 1:numel(modelnames)            
            mfit = ModelBag_get(mbag,modelnames{m});
            if numel(subjs) == 1; mfit = mfit{subjs}; else mfit = mfit(subjs); end
            if fig == 2
                VestBMS_plotBimodalDisparityData(data,task,mfit,Ngen,flags,[],[],[hg_c(m),hg_c(m),hg_c(m)]);
            elseif fig == 3
                VestBMS_plotBimodalPSE(data,task,mfit,Ngen,flags,[],[],hg_c(m));
                set(gca,'Xtick',[-37.5,-25,-10,0,10,25,37.5],'XTickLabel',{'-37.5°','-25°','-10°','0°','10°','25°','37.5°'});
            end
            text(mtxtpos,0.95,modelnames_text{m},'Units','normalized','FontSize',axesfontsize);
            title(''); 
            if m == 1 || m == 3; ylabel(ystring); else ylabel(''); end
            if m < 3; xlabel(''); end
            idx = find(strcmp(g_modelnames,modelnames{m}),1);
            % goftxt = [num2str(mean(g(idx,subjs)),'%.2f') ' ± ' num2str(stderr(g(idx,subjs)),'%.2f')];
            goftxt = [num2str(mean(g(idx,subjs)),'%.2f')];
            text(0.8,0.95,goftxt,'Units','normalized','FontSize',axesfontsize);
            drawnow;
        end
        
        % Plot Bayesian model selection results
        plotfactors(bms{bms_idx}.factor,factornames,hg(2:4),axesfontsize);
        axes(hg(9));
        axis off;
        
        set(gcf,'Position',[1,41,1920,958],'PaperPositionMode','auto');
        saveas(gcf,['fig' num2str(fig) '.svg']);
        saveas(gcf,['fig' num2str(fig) '.png']);   
        
    case 4  % Plot parameters

        % color = [31,120,180; 178,223,138; 166,206,227; 51,160,44]/255;
        color = [178,223,138; 31,120,180; 166,206,227; 0,0,0]/255;
        fontsize = 16;
        plottype = 2;
        
        % Get compatibility data
        C = VestBMS_stats(3,mfits,bms);        
                
        fprintf('\nJoint fits (all datasets):\n');
        mbag_four = {mfits.mbag_uni, mfits.mbag_biml, mfits.mbag_bimu, mfits.mbag_joint};
        modelsummary_four = {mfits.modelsummary_uni, mfits.modelsummary_biml, mfits.modelsummary_bimu, mfits.modelsummary_joint};
        bms_four = {bms{4},bms{1},bms{2},bms{3}};
        params = {'sigma_vest','sigma_vis_low','sigma_vis_med','sigma_vis_high','w_vest','w_vis','lambda','pcommon','kcommon','priorsigma','priorsigmadelta'};
        close all;
        figure;
        for i = 1:7
            subplot(3,4,i);
            [postpdf,pdf,comp.C3(i,:),comp.P3(i,:,:),comp.pxp3(i,:,:)] = ModelPlot_parameterDistribution(mbag_four,params{i},bms_four,[],modelsummary_four,plottype,color);
            cpstring = [num2str(C.pxp3(i,1),'$C_p = %.2f$')];
            if plottype == 1
                text(0.1,1,cpstring,'Units','normalized','FontSize',axesfontsize,'Interpreter','LaTeX');
            else
                xb = xlim;
                xxx = [xb(1) + 0.95*(xb(2) - xb(1)), xb(2), xb(2)];
                text(xb(2),0,cpstring,'FontSize',axesfontsize,'Interpreter','LaTeX','HorizontalAlignment','center');
                plot(xxx,[-3.5 -3.5 0],'k','LineWidth',1);
            end
            legend off;
        end
        mbag_three = {mfits.mbag_biml, mfits.mbag_bimu, mfits.mbag_joint};
        modelsummary_three = {mfits.modelsummary_biml, mfits.modelsummary_bimu, mfits.modelsummary_joint};
        bms_three = {bms{1},bms{2},bms{3}};
        for i = 8:numel(params)
            subplot(3,4,i+1);
            [postpdf,pdf,comp.C3(i,:),comp.P3(i,:,:),comp.pxp3(i,:,:)] = ModelPlot_parameterDistribution(mbag_three,params{i},bms_three,[],modelsummary_three,plottype,color(2:4,:));
            cpstring = [num2str(C.pxp2(i,1),'$C_p = %.2f$')];
            if plottype == 1
                text(0.1,1,cpstring,'Units','normalized','FontSize',axesfontsize,'Interpreter','LaTeX');
            else
                xb = xlim;
                xxx = [xb(1) + 0.95*(xb(2) - xb(1)), xb(2), xb(2)];
                text(xb(2),0,cpstring,'FontSize',axesfontsize,'Interpreter','LaTeX','HorizontalAlignment','center');
                plot(xxx,[-2.5 -2.5 0],'k','LineWidth',1);
            end
            legend off;
        end
        
        for iPanel = [1:7, 9:12]
            subplot(3,4,iPanel);
            params = {'sigma-vest','sigma-vis-low','sigma-vis-med','sigma-vis-high','w-vest','w-vis','lambda','pcommon','kcommon','priorsigma','priorsigmadelta'};
            params_txt = {'$\sigma_\mathrm{vest}$ (deg)','$\sigma_\mathrm{vis-high}$ (deg)','$\sigma_\mathrm{vis-med}$ (deg)','$\sigma_\mathrm{vis-low}$ (deg)', ...
                '$w_\mathrm{vest}$', '$w_\mathrm{vis}$', '$\lambda$','$p_\mathrm{common}$','$k_\mathrm{common}$ (deg)','$\sigma_\mathrm{prior}$ (deg)','$\Delta_\mathrm{prior}$ (deg)'};
            xl = get(gca,'xlabel');
            idx = find(strcmp(xl.String,params),1);
            if isempty(idx) || idx > numel(params_txt); continue; end
            xlabel(params_txt{idx},'Interpreter','LaTeX','FontSize',fontsize);
            if any(iPanel == [1 5 9])
                yl = get(gca,'ylabel');
                ylabel(yl.String,'FontSize',fontsize);
            else
                ylabel('');
            end
        end
        
        subplot(3,4,8);
        for i = 1:4; plot(0,0,'-','LineWidth',3,'Color',color(i,:)); hold on; end
        hl = legend(mfits.mbag_uni.legend, mfits.mbag_biml.legend,mfits.mbag_bimu.legend, mfits.mbag_joint.legend);
        set(hl,'Box','off','Location','West','FontSize',axesfontsize);
        axis off;        
        
        set(gcf,'Position',[1,41,1920,958],'PaperPositionMode','auto');
        saveas(gcf,['fig' num2str(fig) '.svg']);
        saveas(gcf,['fig' num2str(fig) '.png']);   
        
        varargout{1} = comp;        
        
    case 5

        factornames = {'Noise','CI Strategy','Prior'};
        mbag = mfits.mbag_joint;
        modelnames = {'CXD','BPD','BPFs','BPFDs'};
        modelnames_text = {'FIX-X-E','BAY-X-E','BAY/FFU-X-E','BAY/FFU-X-I'};

        Ngen = 20;   % Number of generated datasets per subject
        
        % Full datasets for 11 subjects
        if isempty(subfig); subjs = 1:11; else subjs = subfig; end
        data = data(subjs);
                
        % Plot data
        grid = reshape(1:16,4,4)';
        grid = [1 1 11 11 12 12, 2 2 2 2 2 2; 1 1 11 11 12 12 2 2 2 2 2 2; 1 1 11 11 12 12 2 2 2 2 2 2; 3 3 3 4 4 4 5 5 5 6 6 6;  3 3 3 4 4 4 5 5 5 6 6 6; 3 3 3 4 4 4 5 5 5 6 6 6; 3 3 3 4 4 4 5 5 5 6 6 6; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10];
        hg = plotify(grid,'Margins',[0.1 0.075],'Gutter',[.05, .125],'Labels',{'a','','b','','','','c'},'FontSize',fontsize);        
        
        % Plot model fits
        [g,g_modelnames] = VestBMS_stats([1,3],mfits,bms);
        
        modelsummary = mfits.modelsummary_joint;
        for m = 1:numel(modelnames)
            for iTask = 2:3
                hg_c = hg(2+m+(iTask-2)*4);

                mfit = ModelBag_get(mbag,modelnames{m});
                if numel(subjs) == 1; mfit = mfit{subjs}; else mfit = mfit(subjs); end
                switch iTask
                    case 1
                    
                    
                    case 2
                        VestBMS_plotBimodalDisparityData(data,3,mfit,Ngen,[0 0],[],[],[hg_c,hg_c,hg_c]);
                        ystring = 'Fraction responses ''unity''';
                        legloc = 'NorthWest';   mtxtpos = 0.1;
                        title(modelnames_text{m},'FontSize',fontsize); 
                        if m == 1; ylabel(ystring); else ylabel(''); end
                        % if m < 3; xlabel(''); end
                        idx = find(strcmp(g_modelnames,modelnames{m}),1);
                        goftxt = [num2str(mean(g(idx,subjs)),'%.2f')];
                        text(0.8,1.1,goftxt,'Units','normalized','FontSize',axesfontsize);
                        
                    case 3
                        VestBMS_plotBimodalPSE(data,2,mfit,Ngen,1,[],[],hg_c);
                        set(gca,'Xtick',[-37.5,-25,-10,0,10,25,37.5],'XTickLabel',{'-37.5°','-25°','-10°','0°','10°','25°','37.5°'});
                        mtxtpos = 0.4;
                        ystring = 'Vestibular bias';
                        legloc = 'SouthEast';
                        
                end
                % text(mtxtpos,0.95,modelnames_text{m},'Units','normalized','FontSize',axesfontsize);
                if m == 1; ylabel(ystring); else ylabel(''); end
                % if m < 3; xlabel(''); end
                idx = find(strcmp(g_modelnames,modelnames{m}),1);
                goftxt = [num2str(mean(g(idx,subjs)),'%.2f')];
                % text(0.8,0.95,goftxt,'Units','normalized','FontSize',axesfontsize);
                drawnow;
            end
        end
        
        % Plot Bayesian model selection results
        plotfactors(bms{3}.factor,factornames,hg([1,11,12]),axesfontsize);
        axes(hg(2));
        axis off;
        
        set(gcf,'Position',[1,41,1920,958],'PaperPositionMode','auto');
        saveas(gcf,['fig' num2str(fig) '.svg']);
        saveas(gcf,['fig' num2str(fig) '.png']);           
        
    case 102    % Explicit inference: all subjects
        task = 3;
        flags = [0 0];

        % Model fits
        modelnames = {'BPD','BP-C','CX-C','FF'};
        modelnames_text = {'BAY-X-E','BAY-C-I','FIX-C','SFU'};
        mbag = mfits.mbag_bimu;
        
        % Plot data and fits
        for m = 1:numel(modelnames)
            figure;
            mfit = ModelBag_get(mbag,modelnames{m});
            hg = plotify(3,4,'Margins',[0.1,0.1],'Gutter',[0.04 0.04]);
            for i = 1:Nsubjs
                VestBMS_plotBimodalDisparityData(data(i),task,mfit{i},50,flags,[],[],[hg(i),hg(i),hg(i)]);
                if i == 1; title(modelnames_text{m}); else title(''); end
                if i == 9; xlabel('Stimulus disparity'); else xlabel(''); end
                if i < 9; set(gca,'XTickLabel',''); end
                if i ~= 1 && i ~= 5 && i ~= 9; set(gca,'YTickLabel',''); else ylabel('Proportion response ''unity'''); end
                drawnow;
            end
            axes(hg(Nsubjs+1));
            axis off;
            set(gcf,'Position',[1,41,1920,958],'PaperPositionMode','auto');
            drawnow;
            pause(0.1);
            saveas(gcf,['subjs-unity-' modelnames{m} '.svg']);
            saveas(gcf,['subjs-unity-' modelnames{m} '.png']);
        end
        
        
        
end

end

%--------------------------------------------------------------------------

function plotfactors(factor,factornames,hg,fontsize)
%PLOTFACTORS Plot factor graph
for i = 1:numel(factor)
    axes(hg(i));
    bar(factor{i}.tab,'FaceColor',0.8*[1 1 1],'EdgeColor','none');
    box off;
    set(gca,'TickDir','out','XTickLabel',factor{i}.factornames,'YTick',[0 0.5 1],'FontSize',fontsize);
    axis([0, numel(factor{i}.tab)+1, 0, 1]);
    if i == 1; ylabel('$\tilde{\varphi}$','Interpreter','LaTeX','FontSize',fontsize); end
    xticklabel_rotate([],45);
    if ~isempty(factornames)
        title(factornames{i},'FontSize',fontsize)
    end
end
end
