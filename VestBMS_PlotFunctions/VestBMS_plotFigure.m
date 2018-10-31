function varargout = VestBMS_plotFigure(fig,subfig,mfits,bms,filetype,filesuffix)
%VESTBMS_PLOTFIGURE Plot a specific Figure for the paper or report.
% (Supplementary figures start from 101.)
%
% Usage example:
% VestBMS_plotFigure(fig,data)
% 

fontsize = 24;
axesfontsize = 18;
uppercase = 1;
figpos = [1,41,1920,958];   % Default figure position

Nsubjs = 11;
postflag = 1;           % Plot expected model frequency in BMS plots
Ngen = 20;             % Number of generated datasets per subject (100 for paper)
ticklength = 0.0035;    % Tick length for all panels

if nargin < 2; subfig = []; end

% Load model fit results
if nargin < 3 || isempty(mfits); mfits = load('VestBMS_modelfits'); end
if nargin < 4; bms = []; end

% Run this
% [bms,mfits] = VestBMS_modelComparison('loocv',1,mfits,[],[],[1 1 1 1 0]);

mfits = VestBMS_fixModels(mfits);
mfits.mbag_uni.legend = 'Unisensory discrimination';
mfits.mbag_biml.legend = 'Bisensory implicit CI';
mfits.mbag_bimu.legend = 'Bisensory explicit CI';
mfits.mbag_joint.legend = 'Joint datasets';

% Load datasets
load('VestBMS_data.mat');

plots = VestBMS_defaults('plots');



if uppercase
    panellabels = {'A','B','C','D','E','F','G','H'};
else
    panellabels = {'a','b','c','d','e','f','g','h'};
end

rng(0); % Fix random seed

switch fig
    
    case 1  % Stimuli
        
        fontsize = 36;
        axesfontsize = 28;
        
        bincenters = mfits.mbag_biml.bag{1}.infostruct.bincenters_bim;
        bincenters_vis = bincenters{1};
        bincenters_vest = bincenters{2};

        srange_vis = bincenters_vis(:);
        srange_vest = bincenters_vest(:);
        idx = srange_vis == srange_vest;
        srange_uni = srange_vis(idx);           % C = 1
        srange_vis = srange_vis(~idx);          % C = 2, vis
        srange_vest = srange_vest(~idx);        % C = 2, vest
        col = [0 0 0; 0.6 0.6 0.6];
        % col = [0.6 0.9 0.1; 0.2 0.6 1];
        plot([-50,40],[-50,40],'k--','LineWidth',1); hold on;
        scatter(srange_uni,srange_uni,225,'s','MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:)); hold on;
        scatter(srange_vis,srange_vest,225,'s','MarkerFaceColor',col(2,:),'MarkerEdgeColor',col(2,:));
        xlabel('$s_\mathrm{vis}$','FontSize',fontsize,'Interpreter','LaTeX'); 
        ylabel('$s_\mathrm{vest}$','FontSize',fontsize,'Interpreter','LaTeX');
        axis([-50 50 -50 50]);
        
        ticks = [-45:15:45];
        for i = 1:numel(ticks); ticklabel{i} = [num2str(ticks(i)), '°']; end
        
        hmain = gca;
        set(gca,'XTick',ticks, 'XTickLabel', ticklabel);
        set(gca,'YTick',ticks, 'YTickLabel', ticklabel);
        axis square;
        box off;
        
        set(gca,'TickDir','out','Fontsize',axesfontsize);
        set(gcf,'Color','w');
        set(gcf,'Position',figpos,'PaperPositionMode','auto');
        drawnow;

        % Causal inference panel
        pos = get(gca,'Position');
        %newpos = [pos(1) + pos(3)*0.725, pos(2) + pos(4)*0.1, pos(3)*0.2, pos(4)*0.25];
        newpos = [pos(1) + pos(3)*0.8, pos(2) + pos(4)*0.75, pos(3)*0.2, pos(4)*0.25];
        axes('Position', newpos);
        bar(1,0.2,0.7,'FaceColor',col(1,:),'EdgeColor','none'); hold on;
        bar(2,0.8,0.7,'FaceColor',col(2,:),'EdgeColor','none');
        box off;
        set(gca,'TickDir','out','FontSize',axesfontsize,'XTick',[1 2],'YTick',[0 0.2 0.8]);
        ylabel('Probability','FontSize',fontsize);
        xlabel('$C$','FontSize',fontsize,'Interpreter','LaTeX');
        xlim([0 3]);
        ylim([0 1]);
        
        % Noise panel
        plots = VestBMS_defaults('plots');  % Get default plot info
        col = plots.NoiseColors;            % Colors for different noise levels        
        newpos = [pos(1) + pos(3)*0.8, pos(2) + pos(4)*0.1, pos(3)*0.2, pos(4)*0.25];
        axes('Position', newpos);
        bar(1,1/3,0.7,'FaceColor',col(3,:),'EdgeColor','none'); hold on;
        bar(2,1/3,0.7,'FaceColor',col(2,:),'EdgeColor','none');
        bar(3,1/3,0.7,'FaceColor',col(1,:),'EdgeColor','none');
        box off;
        set(gca,'TickDir','out','FontSize',axesfontsize,'XTick',[1 2 3],'XTickLabel',{'low','med','high'},'YTick',[0 1/3],'YTickLabel',{'0','1/3'});
        ylabel('Probability','FontSize',fontsize);
        xlabel('$c_\mathrm{vis}$','FontSize',fontsize,'Interpreter','LaTeX');
        xticklabel_rotate([],45);
        xlim([0 4]);
        ylim([0 1]);        
        savefigure(['fig' num2str(fig) 'b']);
    
    case 2 % Decision boundaries
        
        if isempty(subfig); subfig = 10; end
        linewidth = 4;
        fontsize = 24;
        axesfontsize = 20;
        
        % Plot midline
        plot([-50,50],[-50,50],'k--','LineWidth',1); hold on;        
        
        % Bayesian model decision rule
        mfit = ModelBag_get(mfits.mbag_joint,'BPD');
        clear functions;
        [~,extras] = ModelWork_like('VestBMS',mfit{subfig});
        for iRel = 1:3        
            xrange_vis = extras.struct{4+iRel}.xrange_vis;
            xrange_vest = extras.struct{4+iRel}.xrange_vest;
            w1 = extras.struct{4+iRel}.w1;
            [~,h(iRel)] = contour(xrange_vis(:), xrange_vest(:), squeeze(w1), [0.5 0.5], 'LineWidth', linewidth, 'Color', plots.NoiseColors(iRel,:)); hold on;            
        end
        
        % Fixed decision rule
        mfit = ModelBag_get(mfits.mbag_joint,'CX');
        % clear functions;
        halfk = mfit{subfig}.mp.fulltheta{4}.kcommon/2;
        h(4) = plot([-100 100] + halfk/sqrt(2)*[1 1], [-100, 100] + halfk/sqrt(2)*[-1 -1], '-.k', 'LineWidth', linewidth); hold on;
        plot([-100 100] + halfk/sqrt(2)*[-1 -1], [-100, 100] + halfk/sqrt(2)*[1 1], '-.k', 'LineWidth', linewidth); hold on;
                
        %[~,extras] = ModelWork_like('VestBMS',mfit{subfig});
        %xrange_vis = extras.struct{4+iRel}.xrange_vis;
        %xrange_vest = extras.struct{4+iRel}.xrange_vest;
        %w1 = extras.struct{4+iRel}.w1;
        %[~,h(iRel)] = contour(xrange_vis(:), xrange_vest(:), squeeze(w1), [0.5 0.5], 'LineWidth', 2, 'Color', plots.NoiseColors(iRel,:)); hold on;            
        
        
        xlabel('$x_\mathrm{vis}$','Interpreter','LaTeX','FontSize',fontsize); 
        ylabel('$x_\mathrm{vest}$','Interpreter','LaTeX','FontSize',fontsize);        
        axis([-50 50 -50 50]);
        ticks = [-45:15:45];
        for i = 1:numel(ticks); ticklabel{i} = [num2str(ticks(i)), '°']; end
        set(gca,'XTick',ticks, 'XTickLabel', ticklabel,'FontSize',axesfontsize);
        set(gca,'YTick',ticks, 'YTickLabel', ticklabel,'TickDir','out');
        axis square; box off;
        set(gcf,'Color','w');
        
        hl = legend(h,'Bayes (High reliability)','Bayes (Medium reliability)','Bayes (Low reliability)','Fixed criterion');
        set(hl,'Location','NorthWest','Box','off','FontSize',fontsize);
        savefigure(['fig' num2str(fig)]);
    
    case {3,4}  % Explicit, implicit inference plot
        
        switch fig
            case 3
                task = 3; stat = 1; bms_idx = 2;
                ystring = 'Fraction responses ''unity''';
                modelnames = {'BPD','BP-C','CX-C','FF'};
                modelnames_text = {'BAY-X-E','BAY-C-I','FIX-C','SFU'};
                mbag = mfits.mbag_bimu;
                legloc = 'NorthWest';   mtxtpos = 0.1;
                factornames = {'Sensory noise','CI strategy','Prior'};
                flags = [0 0];
            case 4
                task = 2; stat = 2; bms_idx = 1;
                % ystring = 'Fraction responses ''right''';
                modelnames = {'FFD','BPD','CXD-C','BPD-C'};
                modelnames_text = {'FFU-X-E','BAY-X-E','FIX-C-E','BAY-C-E'};
                mbag = mfits.mbag_biml;
                mtxtpos = 0.4;
                factornames = {'Sensory noise','CI strategy','Prior'};
                flags(1) = 1;  % Plot bias?
                flags(2) = 1;  % Plot sigma
                if flags(1)
                    ystring = 'Vestibular bias';
                    legloc = 'SouthEast';
                else
                    ystring = 'Vestibular PSE';
                    legloc = 'SouthWest';
                end
        end
        
        % Full datasets for 11 subjects
        if isempty(subfig); subjs = 1:11; else subjs = subfig; end
        data = data(subjs);
                
        % Plot data
        if postflag
            % grid = [1 1 1 1 1 1, 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6;  1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 9 9 9 9 9 9 7 7 7 8 8 8; 9 9 9 9 9 9 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8];
            grid = [repmat([1 1 1 1 1 1, 5 5 5 6 6 6],[10 1]); repmat([1 1 1 1 1 1 7 7 7 8 8 8],[2 1]); repmat([9 9 9 9 9 9 7 7 7 8 8 8],[3 1]); repmat([2 2 3 3 4 4 7 7 7 8 8 8],[5,1])];
        else
            grid = [1 1 1 1 1 1, 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6;  1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 5 5 5 6 6 6; 1 1 1 1 1 1 7 7 7 8 8 8; 9 9 9 9 9 9 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8; 2 2 3 3 4 4 7 7 7 8 8 8];
        end
        hg = plotify(grid,'Margins',[0.1 0.1 0.125 0.05],'Labels',{panellabels{1},panellabels{2},'','',panellabels{3}},'FontSize',fontsize,'Position',figpos);
        hg_a = [hg(1),hg(1),hg(1)];
        if fig == 3
            VestBMS_plotBimodalDisparityData(data,task,[],[],flags,fontsize,axesfontsize,hg_a);
        elseif fig == 4
            VestBMS_plotBimodalPSE(data,task,[],[],flags,fontsize,axesfontsize,hg(1));
        end
        title(''); ylabel(ystring);
        for iCol = 1:3; hcol(iCol) = plot(0,0,'-','Color',plots.NoiseColors(iCol,:),'LineWidth',4); end
        hl = legend(hcol,'High reliability','Medium reliability','Low reliability');
        set(hl,'Box','off','FontSize',axesfontsize,'Location',legloc);
        drawnow;
        
        % Plot inset (Figure 4)
        if fig == 4 && flags(2) == 0
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
            if Ngen == 0; continue; end
            mfit = ModelBag_get(mbag,modelnames{m});
            if numel(subjs) == 1; mfit = mfit{subjs}; else mfit = mfit(subjs); end
            if fig == 3
                VestBMS_plotBimodalDisparityData(data,task,mfit,Ngen,flags,fontsize,axesfontsize,[hg_c(m),hg_c(m),hg_c(m)]);
            elseif fig == 4
                VestBMS_plotBimodalPSE(data,task,mfit,Ngen,flags,fontsize,axesfontsize,hg_c(m));
                % set(gca,'Xtick',[-40,-20,0,20,40],'XTickLabel',{'-40°','-20°','0°','20°','40°'});
                % set(gca,'Xtick',[-37.5,-25,-10,0,10,25,37.5],'XTickLabel',{'-37.5°','-25°','-10°','0°','10°','25°','37.5°'});
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
        plotfactors(bms{bms_idx}.factor,factornames,hg([3 2 4]),axesfontsize,2,[],[],postflag,0);
        axes(hg(9));
        if postflag
            hl = plotfactors(bms{bms_idx}.factor,factornames,hg([9 9 9]),axesfontsize,2,[],0,postflag,1);
            set(hl,'Location','SouthWest');
            title('');
            xlim([-2 -1]);
            text(0.49,0.35,'$\tilde{\varphi}$','Interpreter','LaTeX','FontSize',fontsize,'Units','normalized');
        end
        axis off;
        
        equalticklength(hg,ticklength);
        savefigure(['fig' num2str(fig)]);
        
    case 5  % Plot parameters

        % color = [31,120,180; 178,223,138; 166,206,227; 51,160,44]/255;
        color = [178,223,138; 31,120,180; 166,206,227; 0,0,0]/255;
        fontsize = 24;
        plottype = 2;
        
        % Get compatibility data
        C = VestBMS_stats(3,mfits,bms);        
                
        fprintf('\nJoint fits (all datasets):\n');
        mbag_four = {mfits.mbag_uni, mfits.mbag_biml, mfits.mbag_bimu, mfits.mbag_joint};
        modelsummary_four = {mfits.modelsummary_uni, mfits.modelsummary_biml, mfits.modelsummary_bimu, mfits.modelsummary_joint};
        bms_four = {bms{4},bms{1},bms{2},bms{3}};
        params = {'sigma_vest','sigma_vis_low','sigma_vis_med','sigma_vis_high','w_vest','w_vis','lambda','pcommon','kcommon','priorsigma','priorsigmadelta'};
        logflag = [1,1,1,1,0,0,0,0,1,1,1];
        close all;
        figure;
        for i = 1:7
            subplot(3,4,i);
            [postpdf,pdf,comp.C3(i,:),comp.P3(i,:,:),comp.pxp3(i,:,:),xx] = ModelPlot_parameterDistribution(mbag_four,params{i},bms_four,[],modelsummary_four,plottype,color);
            cpstring = [num2str(C.pxp3(i,1),'$C_p = %.2f$')];
            if plottype == 1
                text(0.1,1,cpstring,'Units','normalized','FontSize',fontsize,'Interpreter','LaTeX');
            else
                xb = xlim;
                xxx = [xb(1) + 0.95*(xb(2) - xb(1)), xb(2), xb(2)];
                text(xb(2),0,cpstring,'FontSize',fontsize,'Interpreter','LaTeX','HorizontalAlignment','center');
                plot(xxx,[-3.5 -3.5 -0.3],'k','LineWidth',1);
            end
            legend off;
            if logflag(i)
                theta = qtrapz(bsxfun(@times, pdf{end}, exp(xx)),2)*diff(xx(1:2));
            else
                theta = qtrapz(bsxfun(@times, pdf{end}, xx),2)*diff(xx(1:2));
            end
            fprintf('%20s = %.2f ± %.2f\n', [params{i} ':'], mean(theta), stderr(theta));            
        end
        mbag_three = {mfits.mbag_biml, mfits.mbag_bimu, mfits.mbag_joint};
        modelsummary_three = {mfits.modelsummary_biml, mfits.modelsummary_bimu, mfits.modelsummary_joint};
        bms_three = {bms{1},bms{2},bms{3}};
        for i = 8:numel(params)
            subplot(3,4,i+1);
            [postpdf,pdf,comp.C3(i,:),comp.P3(i,:,:),comp.pxp3(i,:,:),xx] = ModelPlot_parameterDistribution(mbag_three,params{i},bms_three,[],modelsummary_three,plottype,color(2:4,:));
            cpstring = [num2str(C.pxp2(i,1),'$C_p = %.2f$')];
            if plottype == 1
                text(0.1,1,cpstring,'Units','normalized','FontSize',fontsize,'Interpreter','LaTeX');
            else
                xb = xlim;
                xxx = [xb(1) + 0.95*(xb(2) - xb(1)), xb(2), xb(2)];
                text(xb(2),0,cpstring,'FontSize',fontsize,'Interpreter','LaTeX','HorizontalAlignment','center');
                plot(xxx,[-2.5 -2.5 -0.3],'k','LineWidth',1);
            end
            legend off;
            if logflag(i)
                theta = qtrapz(bsxfun(@times, pdf{end}, exp(xx)),2)*diff(xx(1:2));
            else
                theta = qtrapz(bsxfun(@times, pdf{end}, xx),2)*diff(xx(1:2));
            end
            fprintf('%20s = %.2f ± %.2f\n', [params{i} ':'], mean(theta), stderr(theta));
            if i == 8   % Stats for pcommon
                [h,p,~,tab] = ttest(theta)
                mean(theta)./(tab.sd)
%                 bms_temp = bms{3};                
%                 empirical_flag = cellfun(@(x) any(x == 'd'), lower(bms_temp.modelnames));                
%                 bms_temp.g(:,~empirical_flag) = 0;
%                 [~,pdf,~,~,~,xx] = ModelPlot_parameterDistribution({mfits.mbag_joint},'pcommon',{bms_temp},[],{mfits.modelsummary_joint},0);
%                 theta = qtrapz(bsxfun(@times, pdf{end}, xx),2)*diff(xx(1:2))
            end
        end
                
        hg = [];
        degpanels = [1:4,10:12];
        for iPanel = [1:7, 9:12]
            hg(end+1) = gca;
            subplot(3,4,iPanel);
            % Note that low noise (sigma) corresponds to high reliability, and vice versa
            params = {'sigma-vest','sigma-vis-low','sigma-vis-med','sigma-vis-high','w-vest','w-vis','lambda','pcommon','kcommon','priorsigma','priorsigmadelta'};
            params_txt = {'${\sigma_0}_{\mathrm{vest}}$','$${\sigma_0}_{\mathrm{vis}}(c_\mathrm{high})$','$${\sigma_0}_{\mathrm{vis}}(c_\mathrm{med})$','$${\sigma_0}_{\mathrm{vis}}(c_\mathrm{low})$', ...
                '$w_\mathrm{vest}$', '$w_\mathrm{vis}$', '$\lambda$','$p_\mathrm{c}$','$\kappa_\mathrm{c}$','$\sigma_\mathrm{prior}$','$\Delta_\mathrm{prior}$'};
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
            if any(iPanel == degpanels)
                set(gca,'XTick',log([0.3 1 3 10 30 100]),'XTickLabel',{'0.3°','1°','3°','10°','30°','100°'});
            end
            set(gca,'FontSize',axesfontsize);
        end
        
        subplot(3,4,8);
        for i = 1:4; plot(0,0,'-','LineWidth',3,'Color',color(i,:)); hold on; end
        hl = legend(mfits.mbag_uni.legend, mfits.mbag_biml.legend,mfits.mbag_bimu.legend, mfits.mbag_joint.legend);
        set(hl,'Box','off','Location','West','FontSize',axesfontsize);
        axis off;        
        
        equalticklength(hg,ticklength);
        savefigure(['fig' num2str(fig)]);
        
        varargout{1} = comp;        
        
    case 6  % Joint fits

        factornames = {'Noise','CI Strategy','Prior'};
        mbag = mfits.mbag_joint;
        modelnames = {'CXD','BPD','BPFs','BPFDs'};
        modelnames_text = {'FIX-X-E','BAY-X-E','BAY/FFU-X-E','BAY/FFU-X-I'};
        
        % Full datasets for 11 subjects
        if isempty(subfig); subjs = 1:11; else subjs = subfig; end
        data = data(subjs);
                
        [1 2 2; ...
         1 3 3]
        
        % Plot data
        grid = reshape(1:16,4,4)';
        grid = [1 1 11 11 12 12, 2 2 2 13 13 13; 1 1 11 11 12 12 2 2 2 13 13 13; 1 1 11 11 12 12 2 2 2 13 13 13; 3 3 3 4 4 4 5 5 5 6 6 6;  3 3 3 4 4 4 5 5 5 6 6 6; 3 3 3 4 4 4 5 5 5 6 6 6; 3 3 3 4 4 4 5 5 5 6 6 6; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10; 7 7 7 8 8 8 9 9 9 10 10 10];
        hg = plotify(grid,'Margins',[0.1 0.075],'Gutter',[.05, .125],'Labels',{panellabels{1},'',panellabels{2},'','','',panellabels{3}},'FontSize',fontsize);        
        
        % Plot model fits
        [g,g_modelnames] = VestBMS_stats([1,3],mfits,bms);
        
        modelsummary = mfits.modelsummary_joint;
        for m = 1:numel(modelnames)
            if Ngen == 0; continue; end
            
            for iTask = 2:3
                hg_c = hg(2+m+(iTask-2)*4);

                mfit = ModelBag_get(mbag,modelnames{m});
                if numel(subjs) == 1; mfit = mfit{subjs}; else mfit = mfit(subjs); end
                switch iTask
                    case 1
                    
                    
                    case 2
                        VestBMS_plotBimodalDisparityData(data,3,mfit,Ngen,[0 0],fontsize,axesfontsize,[hg_c,hg_c,hg_c]);
                        ystring = 'Fraction responses ''unity''';
                        legloc = 'NorthWest';   mtxtpos = 0.1;
                        title(modelnames_text{m},'FontSize',axesfontsize); 
                        if m == 1; ylabel(ystring); else ylabel(''); end
                        % if m < 3; xlabel(''); end
                        idx = find(strcmp(g_modelnames,modelnames{m}),1);
                        goftxt = [num2str(mean(g(idx,subjs)),'%.2f')];
                        text(0.8,1.1,goftxt,'Units','normalized','FontSize',axesfontsize);
                        
                    case 3
                        VestBMS_plotBimodalPSE(data,2,mfit,Ngen,1,fontsize,axesfontsize,hg_c);
                        % set(gca,'Xtick',[-37.5,-25,-10,0,10,25,37.5],'XTickLabel',{'-37.5°','-25°','-10°','0°','10°','25°','37.5°'});
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
        plotfactors(bms{3}.factor,factornames,hg([11,1,12]),axesfontsize,2,[],[],postflag,0);
                
        % Plot legend for model comparison
        axes(hg(2));
        if postflag
            hl = plotfactors(bms{3}.factor,factornames,hg([2 2 2]),axesfontsize,2,[],0,postflag,1);
            set(hl,'Location','NorthWest');
            title('');
            xlim([-2 -1]);
            text(0.49,0.35,'$\tilde{\varphi}$','Interpreter','LaTeX','FontSize',fontsize,'Units','normalized');            
        end
        axis off;
        
        % Plot legend for model fits
        axes(hg(13));
        for iCol = 1:3; hcol(iCol) = plot(0,0,'-','Color',plots.NoiseColors(iCol,:),'LineWidth',4); hold on; end
        hl = legend(hcol,'High reliability','Medium reliability','Low reliability');
        set(hl,'Box','off','FontSize',axesfontsize,'Location','SouthWest');
        axis off;
        
        equalticklength(hg,ticklength);        
        savefigure(['fig' num2str(fig)]);
        
    case 7 % Alternative models and fits

        if Ngen > 0
            bms_mrg = VestBMS_modelComparison('marginallike_whmu',1,mfits,[],[],[0 0 1 0 0]);
            bms_baylap = VestBMS_modelComparison('loocv',1,mfits,[],'bayes',[0 0 1 0 1]);
            bms_b1 = VestBMS_modelComparison('loocv',1,mfits,[],[],[0 0 1 0 1]);
            bms_b2 = VestBMS_modelComparison('loocv',1,mfits,[],[],[0 0 1 0 2]);
            close all;
        else
            bms_mrg = bms; bms_baylap = bms; bms_b1 = bms; bms_b2 = bms;
        end
        
        factornames = {'Sensory noise','CI strategy','Prior'};
        grid = [(16:20)', reshape(1:15,[3 5])'];
        for i = 1:size(grid,1); panlab{1+3*(i-1)} = panellabels{i}; end
        panlab = [];
        hg = plotify(grid,'Margins',[0.1 0.1 0.125 0.075],'Labels',panlab,'FontSize',fontsize);
        
        % Plot Bayesian model selection results
        plotfactors(bms{3}.factor,factornames,hg([2 1 3]),axesfontsize,2,1,0,postflag,0);
        plotfactors(bms_mrg{3}.factor,factornames,hg([5 4 6]),axesfontsize,2,0,0,postflag,0);
        plotfactors(bms_baylap{3}.factor,factornames,hg([8 7 9]),axesfontsize,2,0,0,postflag,0);
        plotfactors(bms_b1{3}.factor,factornames,hg([11 10 12]),axesfontsize,2,0,0,postflag,0);
        plotfactors(bms_b2{3}.factor,factornames,hg([14 13 15]),axesfontsize,2,0,1,postflag,0);
        
        axes(hg(16));
        hl = plotfactors(bms{3}.factor,factornames,hg([16 16 16]),axesfontsize,2,[],0,postflag,1);
        set(hl,'Location','NorthWest');
        title('');
        xlim([-2 -1]);        
        
        txt = {'Main','Marginal likelihood','Hyperprior \alpha_0 = 1','Bayesian probability matching (replaced)','Bayesian probability matching (subfactor)'};
        
        for iPanel = 16:20
            axes(hg(iPanel));
            axis off;
            text(0.95,0.5,txt{iPanel-15},'FontSize',fontsize,'Units','Normalized','HorizontalAlignment','right','Interpreter','TeX');
        end
        
        equalticklength(hg,ticklength);        
        savefigure(['fig' num2str(fig)]);

        
    case 101    % Explicit inference, all data
                
        % Full datasets for 11 subjects
        if isempty(subfig); subjs = 1:11; else subjs = subfig; end
        data = data(subjs);
        
        task = 3;
        
        modelnames = {'BPD','CX-C'};
        modelnames_text = {'BAY-X-E','FIX-C'};
        flags = [0 1];
        mbag = mfits.mbag_bimu;
        grid = [[1; 2; 3]*ones(1,10), [8;8;8], [4; 5; 6]*ones(1,10), [7;7;7]];
        hg = plotify(grid,'Margins',[0.075 0.125 0.1 0.1],'Gutter',[0.01 0.1],'Labels',{panellabels{1},'','',panellabels{2}},'FontSize',fontsize);
        
        for m = 1:numel(modelnames)
            idx = (m-1)*3+1;
            hg_c = [hg(idx),hg(idx+1),hg(idx+2)];            
            mfit = ModelBag_get(mbag,modelnames{m});
            if numel(subjs) == 1; mfit = mfit{subjs}; else mfit = mfit(subjs); end
            VestBMS_plotBimodalData(data,task,mfit,Ngen,flags,hg_c,fontsize,axesfontsize);
            axes(hg(idx));
            title(modelnames_text{m},'FontSize',fontsize);
            if m > 1
                set(hg_c(1),'Ylabel',[]);
                set(hg_c(2),'Ylabel',[]);
                set(hg_c(3),'Ylabel',[]);
            end
            set(hg_c(1),'Xlabel',[]);
            set(hg_c(2),'Xlabel',[]);
        end
        
        % Plot legend
        axes(hg(7));
        for iCol = 1:3; hcol(iCol) = plot(0,0,'-','Color',plots.NoiseColors(iCol,:),'LineWidth',4); hold on; end
        hl = legend(hcol,'High reliability','Medium reliability','Low reliability');
        set(hl,'Box','off','FontSize',axesfontsize,'Location','NorthWest');
        axis off;
        
        axes(hg(8)); axis off;
        equalticklength(hg,ticklength);         
        savefigure('figS1');
        
    case 102    % Model recovery
        
        fakedata = load('VestBMS_fakedata_joint.mat');
        mbag_mrec = load('VestBMS_22001.mat');
        modelnames = {'FIX-X-E','BAY-X-E','BAY/FFU-X-I','FIX/FFU-C-I','FIX-X-I','FIX-C-E'};
        
        [recomatrix,err] = ModelWork_modelRecoveryTest(fakedata.data,mbag_mrec.mbag,'aicc',modelnames,fontsize,axesfontsize);
        recomatrix
        
        fprintf('Average percentage of correct recovery: %.1f%%.\n', 100*mean(diag(recomatrix)));
        
        savefigure('figS2',[640 413 600 565]);
        
    case 103    % Explicit inference: all subjects
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
            savefigure(['subjs-unity-' modelnames{m}]);
        end
        
end

end

%--------------------------------------------------------------------------

function hl = plotfactors(factor,factornames,hg,fontsize,ylabidx,titleflag,xticksflag,postflag,legendflag)
%PLOTFACTORS Plot factor graph

if nargin < 6 || isempty(titleflag); titleflag = 1; end
if nargin < 7 || isempty(xticksflag); xticksflag = 1; end
if nargin < 8 || isempty(postflag); postflag = 0; end
if nargin < 9 || isempty(legendflag); legendflag = 0; end

if postflag
    barcol = 0.0*[1 1 1];
    postcol = 0.7*[1 1 1];
else
    barcol = 0.8*[1 1 1];    
end

for i = 1:numel(factor)
    axes(hg(i));
    h(1) = bar(factor{i}.tab,'FaceColor',barcol,'EdgeColor','none');
    
    % Plot exp_r
    if postflag
        hold on;
        for k = 1:size(factor{i}.tab, 2)
            if 0
                h(2) = plot([k-0.4,k+0.4],factor{i}.exp_r(k)*[1 1],'k-','LineWidth',1);
            else
                plot([k k],factor{i}.exp_r(k) + [1 -1]*sqrt(factor{i}.cov_r(k,k)),'-','Color',postcol,'LineWidth',2);
                h(2) = plot([k k],factor{i}.exp_r(k)*[1 1],'ko','LineWidth',2,'MarkerFaceColor',postcol,'MarkerEdgeColor',postcol);
            end
        end
    end
        
    box off;
    set(gca,'TickDir','out','YTick',[0 0.25 0.5 0.75 1],'FontSize',fontsize);
    set(gca,'Ygrid','on');
    axis([0, numel(factor{i}.tab)+1, 0, 1]);
    if i == ylabidx
        if postflag
            ylabel('Probability','FontSize',fontsize);            
        else
            ylabel('$\tilde{\varphi}$','Interpreter','LaTeX','FontSize',fontsize); 
        end
    end
    if ~isempty(factornames) && titleflag
        title(factornames{i},'FontSize',fontsize)
    end
    if xticksflag
        set(gca,'XTickLabel',factor{i}.factornames);
        xticklabel_rotate([],45);
    else
        set(gca,'XTickLabel',[]);
    end
    
    if i == 2 && legendflag
        hl = legend(h,'Protected exceedance probability','Posterior frequency');
        set(hl,'Location','NorthEast','Box','off','FontSize',fontsize);
    end
    
    if postflag
        text(1,0.9,['BOR = ' num2str(factor{i}.bor,'%.2f')],'Fontsize',fontsize,'HorizontalAlignment','right','Units','normalized')
    end
    
end
end

function savefigure(figname,pos)

if nargin < 2 || isempty(pos); pos = [1,41,1920,958]; end

set(gcf,'Position',pos,'PaperPositionMode','auto');
drawnow;
pause(0.1);

saveas(gcf,[figname '.fig']);
try
    saveas(gcf,[figname '.svg']);
    saveas(gcf,[figname '.png']);   
catch
    warning(['Could not save figure ''' figname '''.']);
end

end
