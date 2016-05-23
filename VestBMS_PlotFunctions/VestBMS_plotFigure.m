% CueBMS_PLOTFIGURE Plot a specific Figure for the paper or report.
% (Supplementary figures start from 101.)
%
% Usage example:
% CueBMS_plotFigure(fig, data)
% 
function varargout = VestBMS_plotFigure(fig,subfig,data,mbag,modelsummary,prunemodels,extras,flags,filetype,filesuffix)

if ~exist('data', 'var'); data = []; end
if ~exist('subfig', 'var') || isempty(subfig); subfig = 3; end
if ~exist('mbag', 'var'); mbag = []; end
if ~exist('modelsummary', 'var'); modelsummary = []; end
if ~exist('prunemodels', 'var'); prunemodels = []; end
if ~exist('extras', 'var'); extras = []; end
if ~exist('flags', 'var') || isempty(flags); flags = [0 1]; end
if ~exist('filetype', 'var') || isempty(filetype); filetype = 'pdf'; end
if ~exist('filesuffix', 'var'); filesuffix = []; end

% Keep only Bayesian models
onlybayesian = flags(1);

% Keep the figure on screen at the end
keepfigure = flags(end);

standardunimodalmodels = [ ...
    1 1 1 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  1 MEAN and Gaussian noise/likelihood
    1 1 1 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  2 MEAN and Gaussian noise/single-vis likelihood
    2 2 1 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  3 MEAN and quadratic Gaussian noise/likelihood
    2 2 1 5, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... %  4 MEAN and quadratic Gaussian noise/spatially-constant-likelihood
    2 2 1 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  5 MEAN and quadratic Gaussian noise/single_vis-likelihood
    2 2 1 9, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... %  6 MEAN and quadratic Gaussian noise/all-constant-likelihood
    4 2 1 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  7 MEAN and quadratic Gaussian noise/likelihood (single-w)
    4 2 1 5, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... %  8 MEAN and quadratic Gaussian noise/spatially-constant-likelihood (single-w)
    4 2 1 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... %  9 MEAN and quadratic Gaussian noise/single_vis-likelihood (single-w)
    4 2 1 9, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 10 MEAN and quadratic Gaussian noise/all-constant-likelihood (single-w)
    4 4 1 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 11 MEAN and quadratic Gaussian noise/likelihood (one-w)
    4 4 1 5, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 12 MEAN and quadratic Gaussian noise/spatially-constant-likelihood (one-w)
    4 4 1 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 13 MEAN and quadratic Gaussian noise/single_vis-likelihood (one-w)
    4 4 1 9, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 14 MEAN and quadratic Gaussian noise/all-constant-likelihood (one-w)
    6 8 9 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 15 MEAN and quadratic Gaussian noise/likelihood (free-one-w)
    6 8 9 5, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 16 MEAN and quadratic Gaussian noise/spatially-constant-likelihood (free-one-w)
    6 8 9 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 17 MEAN and quadratic Gaussian noise/single_vis-likelihood (free-one-w)
    6 8 9 9, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 18 MEAN and quadratic Gaussian noise/all-constant-likelihood (free-one-w)
    8 2 1 1, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 15 MEAN and quadratic Gaussian noise/likelihood (1p+1w)
    8 2 1 5, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 16 MEAN and quadratic Gaussian noise/spatially-constant-likelihood (1p+1w)
    8 2 1 7, 1 1 1 1, 1 2 1 1, 1 1 1 1; ... % 17 MEAN and quadratic Gaussian noise/single_vis-likelihood (1p+1w)
    8 2 1 9, 5 1 1 1, 1 2 1 1, 1 1 1 1; ... % 18 MEAN and quadratic Gaussian noise/all-constant-likelihood (1p+1w)
    ];

% Add MEDIAN models
standardunimodalmodels_median = standardunimodalmodels;
standardunimodalmodels_median(:,10) = 3;
standardunimodalmodels = [standardunimodalmodels; standardunimodalmodels_median];

% Standard bimodal models
standardbimodalmodels = [ ...
    2 2 7 1, 1 4 4 2, 1 2 1 3, 3 1 1 1; ... % Bayesian model averaging, vest-vis proportional
    2 2 7 1, 1 4 4 2, 1 2 1 3, 3 2 1 1; ... % Bayesian model selection, vest-vis proportional
    2 2 7 1, 1 4 4 2, 1 2 1 3, 3 3 1 3; ... % Bayesian probability matching, vest-vis proportional        
    2 2 7 1, 1 4 4 2, 1 2 1 3, 3 1 2 1; ... % non-Bayesian 1x, vest-vis proportional
    2 2 7 1, 1 4 4 2, 1 2 1 3, 3 1 3 1; ... % non-Bayesian 3x, vest-vis proportional
    % 2 2 7 1, 1 4 4 2, 1 2 1 3, 3 1 1 2; ... % Bayesian model averaging, vest-vis proportional, free P(C=1) for categorical
    % 2 2 7 1, 1 4 4 2, 1 2 1 3, 3 2 1 2; ... % Bayesian model selection, vest-vis proportional, free P(C=1) for categorical
    % 2 2 7 1, 1 4 4 2, 1 2 1 3, 3 3 1 4; ... % Bayesian probability matching, vest-vis proportional, free P(C=1) for categorical
    ];                

multiplevisflag = (any(modelsummary.cnd == 1) && any(modelsummary.cnd == 2)) || (any(modelsummary.cnd == 5) && any(modelsummary.cnd == 6));

ngen = 30; % Number of generated fake datasets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL PRUNING

if ~isempty(prunemodels)
    if ~iscell(prunemodels); prunemodels = {prunemodels}; end

    for iPrune = 1:length(prunemodels)
    
        switch prunemodels{iPrune}
            case {'const','constant'}     % Prune constant models
                standardunimodalmodels(standardunimodalmodels(:,1) == 1, :) = [];                                
            case {'w1','1w'}           % Prune 1w models  
                standardunimodalmodels(standardunimodalmodels(:,1) == 4 & standardunimodalmodels(:,2) == 4, :) = [];    
            case {'1+1w','2w','w2'}         % Prune 1+1w models  
                standardunimodalmodels(standardunimodalmodels(:,1) == 4 & standardunimodalmodels(:,2) == 2, :) = [];    
            case {'3+1w','4w','w4'}         % Prune 3+1w models  
                standardunimodalmodels(standardunimodalmodels(:,1) == 2 & standardunimodalmodels(:,2) == 2, :) = [];
            case 'w'            % Prune all w models
                standardunimodalmodels(standardunimodalmodels(:,1) == 2 & standardunimodalmodels(:,2) == 2, :) = [];
                standardunimodalmodels(standardunimodalmodels(:,1) == 4 & standardunimodalmodels(:,2) == 2, :) = [];    
                standardunimodalmodels(standardunimodalmodels(:,1) == 4 & standardunimodalmodels(:,2) == 4, :) = [];
            case 'nonbayesian-like'     % Prune non-Bayesian likelihood models
                standardunimodalmodels(standardunimodalmodels(:,4) ~= 1, :) = [];
            case '1x'                   % Remove 1x
                standardbimodalmodels(standardbimodalmodels(:,15) == 2, :) = [];    
            case '3x'                   % Remove 3x
                standardbimodalmodels(standardbimodalmodels(:,15) == 3, :) = [];
            case 'avg'                  % Remove Bayesian model avg
                standardbimodalmodels(standardbimodalmodels(:,14) == 1 & standardbimodalmodels(:,15) == 1, :) = [];
            case 'sel'                  % Remove Bayesian model sel
                standardbimodalmodels(standardbimodalmodels(:,14) == 2, :) = [];
            case 'prm'                  % Remove Bayesian model prm
                standardbimodalmodels(standardbimodalmodels(:,14) == 3, :) = [];
            case 'bayes'                % Remove all Bayesian models
                standardbimodalmodels(standardbimodalmodels(:,14) == 1 & standardbimodalmodels(:,15) == 1, :) = [];
                standardbimodalmodels(standardbimodalmodels(:,14) == 2, :) = [];
                standardbimodalmodels(standardbimodalmodels(:,14) == 3, :) = [];
            case 'nonbayes'             % Remove all non-Bayesian models
                standardbimodalmodels(standardbimodalmodels(:,15) == 2, :) = [];    
                standardbimodalmodels(standardbimodalmodels(:,15) == 3, :) = [];
        end


    end
end

% Keep only MEAN models
standardunimodalmodels(standardunimodalmodels(:,10) ~= 2, :) = [];    


% Using models with lapse?
if any(modelsummary.models(:,13) == 2)
    standardunimodalmodels(:,13) = 2;
    display('Showing unimodal models with lapse.');
end

% Using models with no prior?
if any(modelsummary.models(:,8) == 3)
    standardunimodalmodels(:,8) = 3;
    standardbimodalmodels(:,8) = 3;
    display('Showing models with fixed prior.');
end

% Using models with free one-w?
if ~any(modelsummary.models(:,3) == 9)
    standardunimodalmodels(standardunimodalmodels(:,3) == 9,:) = [];
end

% Using models with sinusoidal eccentricity-dependence of noise?
if any(modelsummary.models(:,1) == 3)
    standardunimodalmodels(standardunimodalmodels(:,1) == 2,1) = 3;
    standardunimodalmodels(standardunimodalmodels(:,2) == 2,2) = 3;
    standardunimodalmodels(standardunimodalmodels(:,1) == 4,1) = 5;
    standardunimodalmodels(standardunimodalmodels(:,2) == 4,2) = 5;
    standardunimodalmodels(standardunimodalmodels(:,1) == 8,1) = 9;
    standardunimodalmodels(standardunimodalmodels(:,2) == 4,2) = 5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFIC DATASETS

% Full bimodal models
fullbimodalmodels = [];
for i = 1:size(standardunimodalmodels, 1)
    addmodels = standardbimodalmodels;
    addmodels(:, 1) = standardunimodalmodels(i, 1);
    addmodels(:, 2) = standardunimodalmodels(i, 2);
    if standardunimodalmodels(i,3) == 9
        addmodels(:, 3) = 9;    % Special free one-w
    end
    addmodels(:, 4) = standardunimodalmodels(i, 4);
    addmodels(:, 5) = standardunimodalmodels(i, 5);
    fullbimodalmodels = [fullbimodalmodels; addmodels];
end

if onlybayesian
    % fullbimodalmodels(fullbimodalmodels(:, 15) ~= 2, :) = [];
    fullbimodalmodels(fullbimodalmodels(:, 15) ~= 1, :) = [];
end

% Font sizes
smallfontsize = 14;
axesfontsize = 16;
fontsize = 18;

colmap = (1 - gray)*0.6 + 0.4;
col = colmap(round((1/4)*(size(colmap, 1)-1)) + 1, :);

% Used for model comparison
switch subfig(1)
    case 1; metric = 'marginallike'; suffix = 'mlike';
    case 2; metric = 'bic'; suffix = 'bic';
    case 3; metric = 'aic'; suffix = 'aic';
end
display(['Using model evidence metric: ' upper(metric) '.']);


switch fig(1)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {2, 3};  % Plot unimodal data (FIG 2: moments; FIG 3: binned)
        % SUBFIG may specify a specific subject
        
        % Command line: load modelfit-unimodal-200x.mat; 
        % VestBMS_plotFigure(2,[3,0],data,mbag_uni,modelsummary_uni)
        
        if isempty(subfig); subfig = [0 0]; end        
        if length(subfig) == 1; subfig(2) = 0; end        
        nid = subfig(2);
        
        if fig == 2
            filename = 'unimod-moments-fit'; 
        elseif fig == 3
            filename = 'unimod-binned-fit'; 
        end
        if nid > 0; filename = [filename '-s' num2str(subfig)]; end
        
        figsize = [1920, 1024];
        
        % Perform model comparison on unimodal data
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,1,standardunimodalmodels,[50, 150],Inf,0,0);        
        bms
        close all;
        
        % Get best models per subject
        for i = 1:size(tab, 1)
            [~,bestmodel] = sort(bms.g(i, :), 2, 'descend');
            mfit{i} = ModelBag_get(mbag,i,models(bestmodel(1), :),{modelsummary.cnd});
        end
        
        if nid > 0; data = data{nid}; mfit = mfit{nid}; end
        
        % Plot unimodal data
        if fig == 2
            gendata = VestBMS_plotUnimodalData(data,mfit,ngen,[],fontsize,axesfontsize);
        elseif fig == 3
            gendata = CueBMS_plotBinnedUnimodalData(data, mfit,ngen,0,fontsize,axesfontsize);            
        end
        
        varargout{1} = bms;
        varargout{2} = gendata;
        varargout{3} = mfit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 4;  % Plot model comparison for unimodal datasets
        
        % Commmand line: load modelfit-unimodal-200x.mat; 
        % CueBMS_plotFigure(4,[],[],[],modelsummary_uni)
                        
        filename = ['unimod-modelcomp-' suffix];
        figsize = [1920, 1024];
        sep = [0., 0.05];        
        alphazero = Inf;
                
        % axesfontsize = 16; fontsize = 20;
        siz = [18 1 7];
        hg = multigraph([ones(1, siz(1)) 2*ones(1, siz(2)) 3*ones(1, siz(3))], sep, [0.05, 0.15]);
                
        axes(hg(2)); axis off; box off;
        
        % Inverted order otherwise problems with the colorbar        
        axes(hg(3));
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,3,standardunimodalmodels,[50, 150],alphazero,1,1);
        
        axes(hg(1));
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,1,standardunimodalmodels,[50, 150],bms,1,1);
        
        h = text(0, 0, '2 log Bayes Factor', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', axesfontsize, 'FontName', 'Arial', 'Interpreter', 'TeX');
%        h = text(0, 0, '$2 \log B$', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', fontsize, 'FontName', 'Arial', 'Interpreter', 'LaTeX');
        set(h, 'Position', [1.14, 0.5, 0], 'Rotation', 270);               

        varargout{1} = models;
        varargout{2} = tab;
        varargout{3} = bms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 5;  % Plot model comparison for bimodal datasets
        
        % Commmand line: load modelfit-bimodal-50x.mat; 
        % [models,tab,bms] = VestBMS_plotFigure(5,[],[],[],modelsummary)
                
        filename = ['bimod-modelcomp-' suffix];     
        figsize = [1920, 1024];
        sep = [0., 0.05];        
        alphazero = Inf;
        % alphazero = 0;
                
        % axesfontsize = 16; fontsize = 20;
        siz = [18 1 7];
        hg = multigraph([ones(1, siz(1)) 2*ones(1, siz(2)) 3*ones(1, siz(3))], sep, [0.05, 0.15]);
                
        axes(hg(2)); axis off; box off;
        
        % Inverted order otherwise problems with the colorbar
        if ~isempty(extras); alphazero = extras; end
        axes(hg(3));
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,3,fullbimodalmodels,[50, 150],alphazero,1,1);
        
        axes(hg(1));
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,1,fullbimodalmodels,[50, 150],bms,1,1);
        
        h = text(0, 0, '2 log Bayes Factor', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', axesfontsize, 'FontName', 'Arial', 'Interpreter', 'TeX');
%        h = text(0, 0, '$2 \log B$', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', fontsize, 'FontName', 'Arial', 'Interpreter', 'LaTeX');
        set(h, 'Position', [1.14, 0.5, 0], 'Rotation', 270);               

        varargout{1} = models;
        varargout{2} = tab;
        varargout{3} = bms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 6;  % Plot model feature comparison for unimodal datasets
        
        % Commmand line: load modelfit-unimodal-200x.mat; 
        % VestBMS_plotFigure(6,[],[],[],modelsummary_uni)
                
        filename = ['unimod-featcomp-' suffix];
        figsize = [1280, 800];
        sep = [0.1, 0.1];        
        alphazero = Inf;
                
        % axesfontsize = 16; fontsize = 20;
        hg = multigraph([1 2 3], sep, [0.1, 0.2]);
                
        % axes(hg(2)); axis off; box off;        
        if ~isempty(extras); alphazero = extras; end

        % External noise model
        axes(hg(1)); modelgroups = [];
        if any(modelsummary.models(:,1) == 3)
            noisemodels = [1 3 5 9];
        else
            noisemodels = [1 2 4 8];            
        end
        modelgroups{1} = standardunimodalmodels(standardunimodalmodels(:, 1) == noisemodels(1) & standardunimodalmodels(:, 12) ~= 4,:);
        modelgroups{2} = standardunimodalmodels(standardunimodalmodels(:, 1) == noisemodels(3) & standardunimodalmodels(:, 2) == noisemodels(3),:);
        modelgroups{3} = standardunimodalmodels(standardunimodalmodels(:, 1) == noisemodels(3) & standardunimodalmodels(:, 2) == noisemodels(2),:);
        modelgroups{4} = standardunimodalmodels(standardunimodalmodels(:, 1) == noisemodels(2),:);
        modelgroups{5} = standardunimodalmodels(standardunimodalmodels(:, 1) == noisemodels(4),:);
        if any(standardunimodalmodels(:, 12) == 4)
            modelgroups{end+1} = standardunimodalmodels(standardunimodalmodels(:, 12) == 4,:);
        end
        
        modelsummary.groupnames = {'Constant', 'Quadratic (1)', 'Quadratic (1+1)', 'Quadratic (3+1)', 'Quadratic (1p+1)', 'Quadratic (m)'};        
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],alphazero,1,1);
        title('External noise model', 'FontSize', fontsize);

        % Internal noise model
        axes(hg(2)); modelgroups = [];
        constantlikelihoods = standardunimodalmodels(:, 1) == 1 | standardunimodalmodels(:, 12) == 4;        
        modelgroups{1} = standardunimodalmodels(standardunimodalmodels(:, 4) == 1 & ~constantlikelihoods,:);
        modelgroups{2} = standardunimodalmodels(standardunimodalmodels(:, 4) == 7 & ~constantlikelihoods,:);
        modelgroups{3} = standardunimodalmodels(standardunimodalmodels(:, 4) == 5,:);
        modelgroups{3} = [modelgroups{3}; standardunimodalmodels(constantlikelihoods & standardunimodalmodels(:, 4) == 1,:)];
        modelgroups{4} = standardunimodalmodels(standardunimodalmodels(:, 4) == 9,:);
        modelgroups{4} = [modelgroups{4}; standardunimodalmodels(constantlikelihoods & standardunimodalmodels(:, 4) == 7,:)];
        modelsummary.groupnames = {'Coher+Ecc', 'Ecc', 'Coher', 'None'};
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],bms,1,1);
        title('Internal noise model', 'FontSize', fontsize);
        ylabel('');
        
        % Decision rule
        axes(hg(3)); modelgroups = [];
        modelgroups{1} = standardunimodalmodels(standardunimodalmodels(:, 10) == 2,:);
        modelgroups{2} = standardunimodalmodels(standardunimodalmodels(:, 10) == 3,:);
        modelsummary.groupnames = {'Mean', 'Median'};
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],bms,1,1);
        title('Decision rule', 'FontSize', fontsize);
        ylabel('');
        
        varargout{1} = models;
        varargout{2} = tab;
        varargout{3} = bms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 7;  % Plot model feature comparison for bimodal datasets
        
        % Commmand line: load modelfit-bimodal-50x.mat; 
        % VestBMS_plotFigure(7,[],[],[],modelsummary)
                
        filename = ['bimod-featcomp-' suffix];
        figsize = [1280, 800];
        sep = [0., 0.1];        
        alphazero = 1;
                
        axesfontsize = 12; fontsize = 16;
        siz = 5;
        hg = multigraph([1*ones(1,siz(1)),0,2*ones(1,siz(1)),0,3*ones(1,siz(1))], sep, [0.1, 0.2]);
                
        % axes(hg(2)); axis off; box off;        
        if ~isempty(extras); alphazero = extras; end

        % External noise model
        axes(hg(1)); modelgroups = [];
        
        normalw = fullbimodalmodels(:,3) == 7;
        modelgroups{1} = fullbimodalmodels(fullbimodalmodels(:, 1) == 1 & fullbimodalmodels(:, 12) ~= 5 & normalw,:);
        modelgroups{2} = fullbimodalmodels(fullbimodalmodels(:, 1) == 4 & fullbimodalmodels(:, 2) == 4 & normalw,:);
        modelgroups{3} = fullbimodalmodels(fullbimodalmodels(:, 1) == 4 & fullbimodalmodels(:, 2) == 2 & normalw,:);
        modelgroups{4} = fullbimodalmodels(fullbimodalmodels(:, 1) == 2 & normalw,:);
        modelsummary.groupnames = {'Constant', 'Quadratic (1)', 'Quadratic (1+1)', 'Quadratic (3+1)'};
        if any(fullbimodalmodels(:,3) == 9)
            modelgroups{5} = fullbimodalmodels(fullbimodalmodels(:, 1) == 6 & fullbimodalmodels(:, 2) == 8 & fullbimodalmodels(:,3) == 9,:);
            modelsummary.groupnames{5} = 'Free quadratic (1)';
        end
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],alphazero,1,1);
        title('External noise model', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);

        % Internal noise model
        axes(hg(2)); modelgroups = [];
        constantlikelihoods = fullbimodalmodels(:, 1) == 1 | fullbimodalmodels(:, 12) == 5;
        modelgroups{1} = fullbimodalmodels(fullbimodalmodels(:, 4) == 1 & ~constantlikelihoods,:);
        modelgroups{2} = fullbimodalmodels(fullbimodalmodels(:, 4) == 7 & ~constantlikelihoods,:);
        modelgroups{3} = fullbimodalmodels(fullbimodalmodels(:, 4) == 5,:);
        modelgroups{3} = [modelgroups{3}; fullbimodalmodels(constantlikelihoods & fullbimodalmodels(:, 4) == 1,:)];
        modelgroups{4} = fullbimodalmodels(fullbimodalmodels(:, 4) == 9,:);
        modelgroups{4} = [modelgroups{4}; fullbimodalmodels(constantlikelihoods & fullbimodalmodels(:, 4) == 7,:)];
        modelsummary.groupnames = {'Coher+Ecc', 'Ecc', 'Coher', 'None'};
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],bms,1,1);
        title('Internal noise model', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);
        ylabel('');
        
        % Causal inference criterion
        axes(hg(3)); modelgroups = []; modelsummary.groupnames = [];
        clumpbayesian = 0;
        if clumpbayesian
            for i = 1:3
                modelgroups{i} = fullbimodalmodels(fullbimodalmodels(:, 15) == i,:); 
            end
            modelsummary.groupnames = {'Bayesian', '1x Criterion', '3x Criterion'};
        else
            groupnames = {'Bayes Avg', 'Bayes Sel', 'Bayes Prm', '1x Criterion', '3x Criterion'};
            for i = 1:3
                temp = fullbimodalmodels(fullbimodalmodels(:, 14) == i & fullbimodalmodels(:, 15) == 1,:);
                if ~isempty(temp)
                    modelgroups{end+1} = fullbimodalmodels(fullbimodalmodels(:, 14) == i & fullbimodalmodels(:, 15) == 1,:);
                    modelsummary.groupnames{end+1} = groupnames{i};
                end
            end
            if any(fullbimodalmodels(:, 15) == 2)
                modelgroups{end+1} = fullbimodalmodels(fullbimodalmodels(:, 15) == 2,:);
                modelsummary.groupnames{end+1} = groupnames{4};
            end
            if any(fullbimodalmodels(:, 15) == 3)
                modelgroups{end+1} = fullbimodalmodels(fullbimodalmodels(:, 15) == 3,:);
                modelsummary.groupnames{end+1} = groupnames{5};
            end
        end
        [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,4,modelgroups,[50, 150],bms,1,1);
        title('Causal inference model', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);
        ylabel('');

        varargout{1} = models;
        varargout{2} = tab;
        varargout{3} = bms;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 8;  % Plot model parameters for unimodal datasets
        
        % Commmand line: load modelfit-unimodal-200x.mat; 
        % CueBMS_plotFigure(8,[],[],mbag_uni,modelsummary_uni,bms)
                
        filename = ['unimod-params-' suffix];
        figsize = [1920, 600];
        sep = [0., 0.1];
                
        axesfontsize = 12; fontsize = 16;
        siz = [3 3 1 1 2 1];
        
        % axes(hg(2)); axis off; box off;
        if ~isempty(extras); bms = extras; else bms = []; end
        if isempty(bms)
            [~,~,bms]=ModelWork_plotModelComparison(modelsummary,metric,3,standardunimodalmodels,[50, 150],Inf,1,1);
            close all;
        end

        % hg = multigraph([1*ones(1,siz(1)),0,2*ones(1,siz(2)),0,3*ones(1,siz(3)),0,4*ones(1,siz(4)),0,5*ones(1,siz(5)),6*ones(1,siz(6))], sep, [0.05, 0.2]);
        hg = multigraph([1*ones(1,siz(1)),0,2*ones(1,siz(2)),0,3*ones(1,siz(3))], sep, [0.05, 0.2]);
        set(gcf,'color','w');
        
        % Parameters
        models = bms.models;
        paramnames = {'sigma_vis_low','sigma_vis_med','sigma_vis_high','sigma_aud',...
            'w_vis_low','w_vis_med','w_vis_high','w_aud', ...
            'priormu', 'priorsigma'};
        for i = 1:length(paramnames)
            paramstat{i} = ModelWork_parameterSummary(mbag,modelsummary,paramnames{i},bms,1,models);
        end
        % Get Weber's fractions only from models with 3+1 and 1+1 w
        for i = 5:8
            paramstat{i} = ModelWork_parameterSummary(mbag,modelsummary,paramnames{i},bms,1,models);
            % paramstat{i} = ModelWork_parameterSummary(mbag,modelsummary,paramnames{i},bms,1,models(models(:,1)==2 | (models(:,1)==4 & models(:,1)==2) ,:));
        end
        % Get Weber's fraction from 1w models
        paramstat{11} = ModelWork_parameterSummary(mbag,modelsummary,'w_vis_low',bms,1,models(models(:,1)==4 & models(:,1)==4,:));
        
        x = 1:length(paramstat);
        for i = 1:length(paramstat)
            % paramstat{i}
            y(i) = paramstat{i}.mean;
            yerr(i) = paramstat{i}.serr;
        end        
        mask = {[1:4],[5:8],[9 10]};
        ylims = [0 20; 0 0.4; -10 40];
        yticks = {0:5:20; 0:0.1:0.4; -10:10:40};
        ystring = {'(deg)', '', '(deg)'};
        param = {{'\sigma_{vis-low}', '\sigma_{vis-med}', '\sigma_{vis-high}', '\sigma_{vest}'}, ...
            {'w_{vis-low}', 'w_{vis-med}', 'w_{vis-high}', 'w_{vest}'}, ...
            {'\mu_{prior}', '\sigma_{prior}'}};
%            {'w_{vis-low}', 'w_{vis-med}', 'w_{vis-high}', 'w_{vest}', 'w_1'}, ...
        
        for i = 1:length(mask)
            axes(hg(i));
            errorbar(1:length(mask{i}), y(mask{i}), yerr(mask{i}), 'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'w');
            axis([0.5 length(mask{i})+0.5, ylims(i, :)]);
            % axis([0.5 length(mask{i})+0.5, ylims(i, :)]);
            box off; 
            set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial', 'TickDir', 'out', 'Ytick', yticks{i});
            set(gca, 'Xtick', 1:length(mask{i}), 'XtickLabel', param{i});
            xticklabel_rotate([],0,[],'Interpreter', 'TeX', 'HorizontalAlignment', 'Center');
            ylabel(ystring{i}, 'FontSize', axesfontsize);
        end
        axes(hg(1));
        title('Base noise parameters', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);        
        axes(hg(2)); hold on;
        title('Weber''s fractions', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);        
        % plot([4.5, 4.5], [0 0.3], 'k:', 'LineWidth', 0.5);
        %axes(hg(4));
        %title('Lapse', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);        
        axes(hg(3)); hold on;
        title('Prior parameters', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.66 1.1]);
        plot([0, length(mask{3})+1], [0 0], 'k:', 'LineWidth', 0.5);
        % axes(hg(6)); set(gca, 'YaxisLocation', 'right');
        

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 9;  % Plot model parameters for bimodal datasets
        
        % Commmand line: load modelfit-bimodal-50x.mat; 
        % VestBMS_plotFigure(9,[],[],mbag,modelsummary,bms)
                
        filename = ['bimod-params-' suffix];
        figsize = [1280, 800];
        sep = [0., 0.1];
                
        axesfontsize = 12; fontsize = 16;
        siz = [3 2 3 2];
                
        if ~isempty(extras); bms = extras; else bms = []; end
        if isempty(bms)
            [~,~,bms]=ModelWork_plotModelComparison(modelsummary,metric,3,fullbimodalmodels,[50, 150],Inf,1,1);
            close all;
        end

        hg = multigraph([1*ones(1,siz(1)),0,2*ones(1,siz(2)),0,3*ones(1,siz(3)),0,4*ones(1,siz(4))], sep, [0.1, 0.2]);
        
        % Parameters
        models = bms.models;
        paramstat{1} = ModelWork_parameterSummary(mbag,modelsummary,'k_vis',bms,1,models);
        paramstat{2} = ModelWork_parameterSummary(mbag,modelsummary,'k_aud',bms,1,models);
        paramstat{3} = ModelWork_parameterSummary(mbag,modelsummary,'pcommon',bms,1,models(models(:,15) == 1,:));
        paramstat{4} = ModelWork_parameterSummary(mbag,modelsummary,'kcommon',bms,1,models(models(:,15) == 2,:));
        paramstat{5} = ModelWork_parameterSummary(mbag,modelsummary,'kcommon',bms,1,models(models(:,15) == 3,:));
        paramstat{6} = ModelWork_parameterSummary(mbag,modelsummary,'kcommon',bms,2,models(models(:,15) == 3,:));
        paramstat{7} = ModelWork_parameterSummary(mbag,modelsummary,'kcommon',bms,3,models(models(:,15) == 3,:));
        paramstat{8} = ModelWork_parameterSummary(mbag,modelsummary,'priormu',bms,1,models);
        paramstat{9} = ModelWork_parameterSummary(mbag,modelsummary,'priorsigma',bms,1,models);
        x = 1:length(paramstat);
        for i = 1:length(paramstat)
            y(i) = paramstat{i}.mean;
            yerr(i) = paramstat{i}.serr;
        end
        
        mask = {[1 2], 3, [4 5 6 7], [8 9]};
        ylims = [0 3; 0 1; 0 30; -20 50];
        yticks = {0:1:3, 0:0.2:1, 0:10:40, -20:10:50};
        ystring = {'Noise multiplier', '(Bayesian) Prior probability', '(Non-Bayesian) Unity criterion (deg)','Prior'};
        param = {{'\gamma_{vis}', '\gamma_{vest}'}, {'p(C=1)'}, {'\Delta_{fix}', '\Delta_{low}', '\Delta_{med}', '\Delta_{high}'}, {'\mu_{prior}','\sigma_{prior}'}};
        for i = 1:length(mask)
            axes(hg(i));
            errorbar(1:length(mask{i}), y(mask{i}), yerr(mask{i}), 'k', 'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'w');
            axis([0.5 length(mask{i})+0.5, ylims(i, :)]);
            % axis([0.5 length(mask{i})+0.5, ylims(i, :)]);
            box off; 
            set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial', 'TickDir', 'out', 'Ytick', yticks{i});
            set(gca, 'Xtick', 1:length(mask{i}), 'XtickLabel', param{i});
            xticklabel_rotate([],0,[],'Interpreter', 'TeX', 'HorizontalAlignment', 'Center');
            ylabel(ystring{i}, 'FontSize', 16);
        end
        axes(hg(2));
        title('Model parameters', 'FontSize', fontsize, 'Units', 'normalized', 'Position', [0.5 1.1]);        
        %axes(hg(end));
        %set(gca, 'YaxisLocation', 'right');
          
        varargout{1} = paramstat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case {10,11,12};  % Plot bimodal data 
        % (FIG 10: vestibular 2D; FIG 11: CAT 2D, FIG 12: vest-CAT disparity)
        % SUBFIG may specify a specific subject
        
        % Command line: load modelfit-bimodal-50x.mat; 
        % VestBMS_plotFigure(10,[3 0],data,mbag,modelsummary,bms)
        
        if isempty(subfig); subfig = [0 0]; end        
        if length(subfig) == 1; subfig(2) = 0; end        
        nid = subfig(2);
        alphazero = 1;
        
        switch fig
            case 10; filename = 'bimod-vest2d-fit'; figsize = [1920, 1024];
            case 11; filename = 'bimod-cat2d-fit'; figsize = [1920, 480];
            case 12; filename = 'bimod-vestcat-fit'; figsize = [1920, 1024];
        end
        if nid > 0; filename = [filename '-s' num2str(nid)]; end
        
        % Perform model comparison on bimodal data
        fullmbag = [];
        bms = [];
        if ~isempty(extras)
            if isfield(extras,'g')
                bms = extras;
                models = bms.models;
                tab = bms.g;
            elseif isfield(extras,'bag')
                fullmbag = extras;
            end
            clear extras;
        end
        if isempty(bms)
            [models,tab,bms]=ModelWork_plotModelComparison(modelsummary,metric,1,fullbimodalmodels,[50, 150],alphazero,0,0);
            bms
            close all;
        end
        
        if isfield(mbag, 'prefix')
            % Get best models per subject
            for i = 1:size(tab, 1)
                [~,bestmodel] = sort(bms.g(i, :), 2, 'descend');
                mfit{i} = ModelBag_get(mbag,i,models(bestmodel(1), :),{modelsummary.cnd});
                
                if ~isempty(fullmbag) 
                    % Take data from mbag with full data matrix
                    temp = ModelBag_get(fullmbag,i,models(bestmodel(1), :),{modelsummary.cnd});
                    mfit{i}.X = temp.X;
                end                
            end

            if nid > 0; data = data{nid}; mfit = mfit{nid}; end
        else
            mfit = mbag;
        end
        
        % Plot unimodal data
        switch fig 
            case 10
                [fig,gendata] = VestBMS_plotBimodalData(data,2,mfit,ngen,0,fontsize,axesfontsize);
            case 11
                [fig,gendata] = VestBMS_plotBimodalData(data,3,mfit,ngen,0,fontsize,axesfontsize);
            case 12
                [fig,gendata] = VestBMS_plotBimodalDisparityData(data,[2 3],mfit,ngen,0,fontsize,axesfontsize);
        end
        
        varargout{1} = bms;
        varargout{2} = mfit;
        varargout{3} = fig;
        varargout{4} = gendata;
end

% Save figure and go home
set(gcf,'PaperPositionMode','auto','color','w');
set(gcf, 'Visible', 'on', 'Position', [1, 1, figsize]);
export_fig(filetype, 'm2', 'transparent', [filename filesuffix '.' filetype], gcf);
if ~keepfigure; close all; end

end