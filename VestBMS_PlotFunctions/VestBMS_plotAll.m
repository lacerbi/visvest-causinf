% VestBMS_PLOTALL Plot all figures.
%
function VestBMS_plotAll(plotUnimodal,plotBimodal)

if ~exist('plotUnimodal', 'var') || isempty(plotUnimodal); plotUnimodal = 0; end
if ~exist('plotBimodal', 'var') || isempty(plotBimodal); plotBimodal = 1; end


uni_file = 'vestbms-mfit-uni-500x.mat';
bim_files = {'vestbms-mfit-bim-50x.mat','vestbms-mfit-bim-estonly-50x.mat','vestbms-mfit-bim-catonly-50x.mat'};
filetype = 'png';

flags = [1 1 0; 0 1 0];

load vestdata.mat;

if plotUnimodal
    display('Plotting unimodal data');
    % Load unimodal data
    unistuff = load(uni_file);

    for iFlag = 1:size(flags,1)
        % Plot unimodal data
        CueBMS_plotFigure(2,0,data,unistuff.mbag_uni,unistuff.modelsummary_uni,[],flags(iFlag,:),filetype);
        CueBMS_plotFigure(3,0,data,unistuff.mbag_uni,unistuff.modelsummary_uni,[],flags(iFlag,:),filetype);

        % Model comparison and analysis for marginal likelihood, BIC and AIC
        for iMethod = 1:3
            [~,~,bms] = CueBMS_plotFigure(4,iMethod,[],unistuff.mbag_uni,unistuff.modelsummary_uni,[],flags(iFlag,:),filetype);
            CueBMS_plotFigure(6,iMethod,[],unistuff.mbag_uni,unistuff.modelsummary_uni,bms,flags(iFlag,:),filetype);
            CueBMS_plotFigure(8,iMethod,[],unistuff.mbag_uni,unistuff.modelsummary_uni,bms,flags(iFlag,:),filetype);
        end
    end
    clear unistuff;
end

if plotBimodal
    display('Plotting bimodal data');
    % Load bimodal data
    for ii = 1:length(bim_files)
        bimstuff{ii} = load(bim_files{ii});
    end
    
    methodname = {'aic','bic'};
    prunemodels = {{'3x','bayes'},{'bayes'},{'1x','3x'},[],{'sel','prm','3x'}};
    modelsuffix = {'-only1x','-nonbayesian','-bayesian','-all','-avg_1x'};
    fitsuffix = {'','-estonly','-catonly'};
    methodnumber = [3 2];

    % Model comparison and analysis
    for iMethod = 1:length(methodnumber)
        subfig = methodnumber(iMethod);
        for iModels = 1:length(modelsuffix)
            for iFit = 1:length(fitsuffix)
                suffix = [fitsuffix{iFit} modelsuffix{iModels}];            
                VestBMS_plotFigure(7,subfig,[],[],bimstuff{iFit}.modelsummary, ...
                    prunemodels{iModels},[],[],filetype,suffix);
                close all;
                VestBMS_plotFigure(9,subfig,[],bimstuff{iFit}.mbag,...
                    bimstuff{iFit}.modelsummary,prunemodels{iModels},[],[],filetype,suffix);
                close all;
                VestBMS_plotFigure(12,subfig,data,bimstuff{iFit}.mbag,...
                    bimstuff{iFit}.modelsummary,prunemodels{iModels},bimstuff{1}.mbag,...
                    [],filetype,['-' methodname{iMethod} suffix]);
                close all;
            end
        end
    end
end
end