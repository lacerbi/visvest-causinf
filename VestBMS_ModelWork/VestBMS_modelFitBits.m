

function varargout = VestBMS_modelFitBits(command, varargin)

switch lower(command)

    % function [datalikeFun, infostruct] = ModelDependentDefine(dataone)
    % MODELDEPENDENTDEFINE Define model-dependent variables.
    case 'define'
        dataone = varargin{1};
        model = varargin{2};
        options = varargin{3};

        % Define data function
        datalikeFun = @VestBMS_datalike;        
        varargout{1} = datalikeFun;
        
        if ~isempty(options)             
            cnd = options.cnd;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Specific data for model-parameters structure
            infostruct.experimentName = options.experimentName;
            if isempty(infostruct.experimentName)
                error('Expriment name not specified.');
            end
            
            % 'Temperature' during slice sampling (multiplies window step)
            infostruct.samplingtemperature = options.samplingtemperature;

            % BINNED STIMULI
            
            infostruct.bincenters = options.bincenters(:);
            
            % Bin weights for stimuli (non-Trevor datasets)
            infostruct.binweights = options.binweights;
            
            % Response bins
            infostruct.respbincenters = options.respbincenters;
            if ~isempty(infostruct.bincenters) && isempty(infostruct.respbincenters)
                infostruct.respbincenters = infostruct.bincenters;
            end
            
            if ~isempty(infostruct.bincenters)
                if isempty(options.bincenters_bim)
                    infostruct.bincenters_bim{1} = repmat(infostruct.bincenters, [length(infostruct.bincenters), 1]);
                    infostruct.bincenters_bim{2} = repmat(infostruct.bincenters', [length(infostruct.bincenters), 1]);
                    infostruct.bincenters_bim{2} = infostruct.bincenters_bim{2}(:);
                    infostruct.bincenters_bim{3} = infostruct.respbincenters;
                else
                    infostruct.bincenters_bim = options.bincenters_bim;                    
                end
            else
                infostruct.bincenters_bim = [];
            end
            
            % Maximum stimulus range
            if model(9) == 1    % Uncorrelated stimuli
                infostruct.MAXRNG = 90;
            else                % Correlated stimuli
                infostruct.MAXRNG = 60;
                infostruct.MAXDELTA = 60;
            end
            
            if options.loadinitfromconst && numel(options.dataid) > 1
                model_const = model;
                model_const([1 2]) = 1; % Constant noise model
                
                if options.dataid(2) == 4
                    temp = load('VestBMS_starting_points_unity');
                elseif options.dataid(2) == 8
                    temp = load('VestBMS_starting_points_localization');
                else
                    temp = [];
                end
                if ~isempty(temp)
                    mfit = ModelBag_get(temp.mbag,options.dataid,model_const,cnd);                    
                end
                if ~isempty(mfit) && ~isempty(mfit.mp)
                    fprintf('Loading initial point from constant noise model.\n');
                    Nrep = 1e3;
                    maptheta_const = mfit.maptheta;
                    params_const = mfit.mp.params;
                    infostruct.cnd = cnd;
                    mp = VestBMS_setupModel([],[],model,infostruct);
                    params = mp.params;
                    x0 = NaN(Nrep,numel(params));
                    for i = 1:numel(params_const)
                        idx = find(strcmp(params_const{i},params),1);
                        % Small perturbation
                        x0(:,idx) = maptheta_const(i) + 0.1*mp.bounds.SCALE(idx)*randn(Nrep,1);
                        if isfield(options,'replica') && options.replica == 1; x0(1,idx) = maptheta_const(i); end
                    end
                    % Set eccentricity parameters to small values, others to reasonable values
                    for i = find(isnan(x0(1,:)))
                        if params{i}(1) == 'w'
                            w = 0.2*linspace(1/(Nrep),1,Nrep);
                            x0(:,i) = w(randperm(Nrep));
                            if isfield(options,'replica') && options.replica == 1; x0(1,i) = 0.01; end
                        else
                            x0(:,i) = 0.5*(mp.bounds.RLB(i) + mp.bounds.RUB(i));
                        end
                    end
                    x0 = bsxfun(@min, bsxfun(@max, x0, mp.bounds.LB + 1e-4), mp.bounds.UB - 1e-4); 
                    options.startx = x0;
                end
                    
            end
            
            if options.loadinitfromdisc && numel(options.dataid) > 1
                model_disc = model;
                % Set discrete prior model
                model_disc(9) = model_disc(9) + 2;
                if model_disc(9) < 4 || model_disc(9) > 5; model_disc(9) = 5; end
                if options.dataid(2) == 4
                    temp = load('VestBMS_starting_points_unity');
                elseif options.dataid(2) == 8
                    temp = load('VestBMS_starting_points_localization');
                else
                    temp = [];
                end
                if ~isempty(temp)
                    mfit = ModelBag_get(temp.mbag,options.dataid,model_disc,cnd);                    
                end
                if ~isempty(mfit) && ~isempty(mfit.mp)
                    fprintf('Loading initial point from discrete-prior model.\n');
                    infostruct.cnd = cnd;
                    mp = VestBMS_setupModel([],[],model,infostruct);
                    x0 = mfit.maptheta;
                    x0 = bsxfun(@min, bsxfun(@max, x0, mp.bounds.LB + 1e-4), mp.bounds.UB - 1e-4);
                    options.startx = x0;
                end                    
            end
            
            
            
            % Parameter structure from unimodal data
            if ~isempty(options.unifulltheta)
                infostruct.unifulltheta = options.unifulltheta;
            elseif isfield(dataone, 'unifulltheta')
                infostruct.unifulltheta = dataone.unifulltheta{options.unifullthetanumber};
            else
                infostruct.unifulltheta = [];
            end
            
            % Optimization with fake data starting near the true value
            if (options.fakedatastartx)
                standardbimodalmodels = [ ...
                    2 2 6 1, 1 4 4 6, 1 2 1 3, 3 1 1 1; ... % Bayesian model averaging, audio proportional
                    2 2 6 1, 1 4 4 6, 1 2 1 3, 3 2 1 1; ... % Bayesian model selection, audio proportional
                    2 2 6 1, 1 4 4 6, 1 2 1 3, 3 3 1 3; ... % Bayesian probability matching, audio proportional        
                    2 2 6 1, 1 4 4 6, 1 2 1 3, 3 1 2 1; ... % non-Bayesian 1x, audio proportional
                    2 2 6 1, 1 4 4 6, 1 2 1 3, 3 1 3 1; ... % non-Bayesian 3x, audio proportional
                    ];

                unimodal2bimodalmodels = [ ...
                    4 4 1 1, 1 1 2 5, 1 2 1 2, 2 1 1 1; ... % 11 MEAN and quadratic Gaussian noise/likelihood
                    4 4 1 5, 5 1 2 5, 1 2 1 2, 2 1 1 1; ... % 12 MEAN and quadratic Gaussian noise/spatially-constant-likelihood
                    4 4 1 7, 1 1 2 5, 1 2 1 2, 2 1 1 1; ... % 13 MEAN and quadratic Gaussian noise/single_vis-likelihood
                    4 4 1 9, 5 1 2 5, 1 2 1 2, 2 1 1 1; ... % 14 MEAN and quadratic Gaussian noise/all-constant-likelihood
                    ];

                models = [];
                for i = 1:size(unimodal2bimodalmodels, 1)
                    addmodels = standardbimodalmodels;
                    addmodels(:, 1) = unimodal2bimodalmodels(i, 1);
                    addmodels(:, 2) = unimodal2bimodalmodels(i, 2);
                    addmodels(:, 3) = 7; % Audio-video independently proportional from data
                    addmodels(:, 4) = unimodal2bimodalmodels(i, 4);
                    addmodels(:, 5) = unimodal2bimodalmodels(i, 5);
                    addmodels(:, 12) = unimodal2bimodalmodels(i, 12) + 1;
                    models = [models; addmodels];
                end                
                
                modelnumber = find(all(bsxfun(@eq, models, model), 2));
                startx = dataone.startx{modelnumber};
                options = setoptions(options,'startx',[NaN(size(startx)); startx]);                
            end

            % Use some analytical approximation in the bimodal trials 
            infostruct.approxflag = options.approxflag;

            % Convert model number to model vector
            % if isscalar(model); model = CueBMS_getmodelnumber(model); end
            
            varargout{2} = infostruct;
            varargout{3} = options;
        end
                
    % function modelstruct = ModelDependentInit(dataone, modelstruct, options)
    %  Model-dependent preprocessing of data structure
    case 'preprocessdata'
        
        modelstruct = varargin{2};
        options = varargin{3};
        infostruct = varargin{4};
                
        modelstruct.X = varargin{1}.X;                        
        modelstruct.nData = 0;
        
    %----------------------------------------------------------------------
    % DATASET FLAGS: Datasets keep/remove trials
    % Flag 1: Keep 1st session (0), Skip 1st session (1)
    % Flag 2: Include visual estimation trials (0), Exclude them (1)
    % Flag 3: Include audio estimation trials (0), Exclude them (1)
    % Flag 4: Include categorization trials (0), Exclude them (1)
    % Flag 5: Include zero disparities (0), Exclude them (1)
    % Flag 6: Include large disparities (0), Exclude them (1)
    % Flag 7: Include unisensory trials (0), Exclude them (1)
    %----------------------------------------------------------------------
        
        options.dataid = [options.dataid 0];

        % Model-dependent data processing
        dataflags = num2flags(options.dataid(2),'l',7);
        
        % Remove trials with NaNs
        for iNoise = 1:4
            modelstruct.X.unimodal{iNoise}(any(isnan(modelstruct.X.unimodal{iNoise}),2), :) = [];
        end
        for iNoise = 1:3
            for iType = 1:3
                modelstruct.X.bimodal{iNoise}{iType}(any(isnan(modelstruct.X.bimodal{iNoise}{iType}),2), :) = [];
            end
            modelstruct.X.bimodalall{iNoise}(any(isnan(modelstruct.X.bimodalall{iNoise}),2), :) = [];
        end
        
        % Fix categorical trials with response 0 to 2
        for iNoise = 1:3
            modelstruct.X.bimodal{iNoise}{3}(modelstruct.X.bimodal{iNoise}{3}(:,end) == 0, end) = 2;
            f = modelstruct.X.bimodalall{iNoise}(:,2) == 3 & modelstruct.X.bimodalall{iNoise}(:,end) == 0;
            modelstruct.X.bimodalall{iNoise}(f,end) = 2;
        end
                        
        % If specified, keep only one half of the trials
        %trialMax = modelstruct.X.all(end, 1);
        %switch options.keeponlyhalf
        %    case 0 % Keep all the data
        %        f = 1:trialMax;
        %    case 1 % Keep only first half of all data
        %        f = 1:floor(trialMax/2);
        %    case 2 % Keep only second half of all data
        %        f = floor(trialMax/2)+1:trialMax;
        %end
        
        trialMax = modelstruct.X.all(end, 1);
        f = 1:trialMax;
                
        if dataflags(1)     % Trim the first n trials
            f(f <= options.trimfirstntrials) = [];
        end
        
        % Keep only some sessions
        if isfield(options,'sessions') && ~isempty(options.sessions)
            f = intersect(f,find(any(bsxfun(@eq, modelstruct.X.all(:,2), options.sessions), 2)));
        end

        % Remove flagged trials
        fun = @(X) all(~bsxfun(@eq, X(:,1), f),2);
        modelstruct.X = removetrials(modelstruct.X,fun);

        % Exclude trials of a given type 
        % (visual estimation, audio estimation, categorical)
        for iFlag = 2:4     
            if dataflags(iFlag)
                iType = iFlag - 1;
                if iType == 1
                    for iNoise = 1:3; modelstruct.X.unimodal{iNoise} = []; end
                elseif iType == 2
                    modelstruct.X.unimodal{4} = [];
                end                    
                for iNoise = 1:3
                    modelstruct.X.bimodal{iNoise}{iType} = [];
                    modelstruct.X.bimodalall{iNoise}(modelstruct.X.bimodalall{iNoise}(:, 2) == iType, :) = [];                
                end                    
            end
        end
        
        % Exclude unisensory trials
        if dataflags(7)
            for iNoise = 1:4; modelstruct.X.unimodal{iNoise} = []; end
        end

        % Remove trials with SMALL disparity (5 deg or less)
        if dataflags(5)
            DISPARITY = 5;
            % Consider the case of bimodal and bimdisp
            fun = @(X) logical((size(X,2) == 5) .* (abs(X(:, 3) - X(:, 4)) <= DISPARITY) ...
                + (size(X,2) == 4) .* (abs(X(:, 3)) <= DISPARITY));
            modelstruct.X = removetrials(modelstruct.X,fun,'bimodal');
        end

        % Remove trials with LARGE disparity (10 deg or more)
        if dataflags(6)
            DISPARITY = 10;
            % Consider the case of bimodal and bimdisp
            fun = @(X) logical((size(X,2) == 5) .* (abs(X(:, 3) - X(:, 4)) >= DISPARITY) ...
                + (size(X,2) == 4) .* (abs(X(:, 3)) >= DISPARITY));
            modelstruct.X = removetrials(modelstruct.X,fun,'bimodal');
        end

        % Keep or discard noise levels depending on condition
        for iNoise = 1:3
            if all(options.cnd ~= (iNoise+4))
                modelstruct.X.bimodalall{iNoise} = [];
                for iType = 1:3; modelstruct.X.bimodal{iNoise}{iType} = []; end
            end
        end
            
        % Count unimodal trials
        for iNoise = 1:4
            if any(options.cnd == iNoise)
                modelstruct.nData = modelstruct.nData + size(modelstruct.X.unimodal{iNoise}, 1);                
            end
        end
        
        % Count bimodal trials
        for iNoise = 1:3
            modelstruct.nData = modelstruct.nData + size(modelstruct.X.bimodalall{iNoise}, 1);
        end
        
        % Prepare stimulus/response bin counts for non-Trevor datasets
        if ~isempty(infostruct.bincenters)
            respbincenters = infostruct.respbincenters;
            
            % Unimodal data
            DD = modelstruct.X.unimodal;
            x = [];
            for iicnd = 1:4
                if ~isempty(DD{iicnd}); x = [x; DD{iicnd}(:,2)]; end
            end                
            infostruct.bincenters_uni = unique(x);
            
            for iicnd = 1:length(DD)
                xmat = zeros(numel(infostruct.bincenters_uni),numel(respbincenters));
                for iBin = 1:numel(infostruct.bincenters_uni)
                    for jBin = 1:length(respbincenters)
                        if ~isempty(DD{iicnd})
                            xmat(iBin,jBin) = sum(DD{iicnd}(:, 2) == infostruct.bincenters_uni(iBin) & DD{iicnd}(:, 3) == respbincenters(jBin));
                        end
                    end
                end
                modelstruct.X.unibins{iicnd} = xmat;
            end
            
            % Bimodal data
            bincenters_vis = infostruct.bincenters_bim{1};
            bincenters_vest = infostruct.bincenters_bim{2};
            binmeshLength = size(bincenters_vis,1);
            DD = modelstruct.X.bimodal;
            
            for iNoise = 1:length(DD)
                for iType = 1:3
                    if iType == 3
                        thisrespbincenters = [1 2]; % Unity
                    else
                        thisrespbincenters = respbincenters; % Estimation
                    end                    
                    if isempty(DD{iNoise}{iType})
                        xmat = [];
                    else
                        xmat = zeros(binmeshLength,length(thisrespbincenters));
                        for iBin = 1:binmeshLength
                            for jBin = 1:length(thisrespbincenters)
                                xmat(iBin,jBin) = sum(DD{iNoise}{iType}(:, 4) == bincenters_vis(iBin) & ...
                                DD{iNoise}{iType}(:, 3) == bincenters_vest(iBin) & ...
                                    DD{iNoise}{iType}(:, end) == thisrespbincenters(jBin));
                            end
                        end
                    end
                    modelstruct.X.bimbins{iNoise}{iType} = xmat;
                end
            end
            
        end
        
        varargout{1} = modelstruct;
        varargout{2} = infostruct;
        
    % function cost = ModelDependentCost(dataone, model)
    % MODELCOST Model-dependent cost.
    case 'cost'
        
        model = varargin{1};
        cnd = varargin{2};
               
        % Compute approximate model computational cost
        % priornmix = [1 1 2 2];
        % lapsecost = [1 1.25 1.5 1.25 1.5 1.25 1.25];
        cost = ones(1, length(cnd));
        for j = 1:length(cnd)
           % ...        
        end   
        
        varargout{1} = sum(cost);
        
end

end

%--------------------------------------------------------------------------
function X = removetrials(X,fun,type)
%REMOVETRIALS Remove trials from model data X according to rule FUN
%   TYPE can be 'all', 'unimodal', 'bimodal'.

if nargin<3; type = 'all'; end

if strcmpi(type(1),'a')
    X.all(fun(X.all), :) = [];
end
   
if strcmpi(type(1),'a') || strcmpi(type(1),'u')
    for iNoise = 1:4
        if ~isempty(X.unimodal{iNoise})
            X.unimodal{iNoise}(fun(X.unimodal{iNoise}), :) = [];
        end
    end
end

if strcmpi(type(1),'a') || strcmpi(type(1),'b')
    for iNoise = 1:3
        for iType = 1:3
            if ~isempty(X.bimodal{iNoise}{iType})
                X.bimodal{iNoise}{iType}(fun(X.bimodal{iNoise}{iType}), :) = [];
            end
            if ~isempty(X.bimdisp{iNoise}{iType})        
                X.bimdisp{iNoise}{iType}(fun(X.bimdisp{iNoise}{iType}), :) = [];
            end
        end
        if ~isempty(X.bimodalall{iNoise})
            X.bimodalall{iNoise}(fun(X.bimodalall{iNoise}), :) = [];
        end
    end
end

end
