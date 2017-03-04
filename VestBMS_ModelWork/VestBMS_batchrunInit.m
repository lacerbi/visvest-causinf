function [options,models,dataids,cnd] = VestBMS_batchrunInit(data,type,options)
% VESTBMS_BATCHRUNINIT initialize variables for batch run.

if nargin < 3; options = []; end

debug = 0;

% Get additional type parameters
type = type(1);

% Trials are binned for computation of the log likelihood (obsolete!)
% options = setoptions(options,'binnedloglik',1,1);

% Continue previous sampling if exists
options = setoptions(options,'loadstartx',1,1);

% Slice sampling window multiplier
options = setoptions(options,'samplingtemperature',1,1);

% By default do not compute Hessian unless optimizing (see below)
options = setoptions(options,'hessianflag',0,1);

options = setoptions(options,'optimizationmethod','bps',1);
% options = setoptions(options,'optimizationmethod','fmincon',1);

% Get host name
[~,hostname] = system('hostname');
options = setoptions(options,'hostname',strtrim(hostname),1);

% Number of datasets
if isfield(data,'data'); nDatasets = length(data.data);
else nDatasets = length(data); end

% Subjects mask
DATAIDS = [(1:nDatasets)',zeros(nDatasets,1)];

% Remove monkey data
DATAIDS(any(bsxfun(@eq, DATAIDS(:,1), [12 13 14]),2),:) = [];

% Default number of samples for unimodal/bimodal trials
NSAMPLES = [1e4,5e3];
NSAMPLES = [0 0];

% Default optimization steps before starting sampling
MAXFUNEVALS = 1500;

% Optimization steps when optimizing only
NITER_OPTIMIZATION = 1500;

% Number of restarts for optimization
nOptimizationRestarts = 50;

if debug
    nOptimizationRestarts = 10;
    NSAMPLES = [10,10];
    MAXFUNEVALS = 10;
end

options = setoptions(options,'nstarts',nOptimizationRestarts,1);
options = setoptions(options,'nsobol',1e4,1);
options = setoptions(options,'optfevals',MAXFUNEVALS,1);
if isempty(options.optfevals); options.optfevals = MAXFUNEVALS; end

dataids = DATAIDS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT MODELS

models = [];
standardunimodalmodels = [ ...
    5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 1 1 0; ... % Base
    ];

%unimodal2bimodalmodels = [ ...
%    5 3 1 1, 1 1 1 3, 1 1 2 1, 1 1 2 1 0; ... %  1 MEAN and Gaussian noise/likelihood
%    ];



% Best unimodal models
bestunimodalmodels = [ ...
    4 2 1 5, 5 1 2 5, 1 2 1 2, 2 1 1 1 0; ... %  9 MEAN and quadratic Gaussian noise/spatially-constant-likelihood (single w)
    ];

% Standard bimodal models
standardbimodalmodels = [ ...
    5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 2 1 0; ... % Generalized Bayesian
    5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 4 1 0; ... % Soft fixed criterion
    5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 5 1 0; ... % Forced fusion
    ];


switch type
    case 0; % DEBUG    
        
        [options,models,groupcnd] = VestBMS(options,2,0);
        options.jobname = 'vest_debug';
        dataids = [1 0];
        
        options = setoptions(options,'nsamples',100,1);
        options = setoptions(options,'nstoredsamples',15,1);
        options = setoptions(options,'optfevals',20,1);
        options = setoptions(options,'nstarts',1,1);
        models = models(1,:);
        
%--------------------------------------------------------------------------        
% BISENSORY ESTIMATION DATA FITS
    
   case {1} % Bisensory standard models
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_bim';
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
    case 2  % Monkey only
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_bim_monkey';
       dataids = [12 8; 13 8; 14 8];
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case {11} % Bisensory standard models without beta
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_2bim';
       models(:,11) = 3;    % Probability matching
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case {21} % Bisensory standard models without beta and with lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_lbim';
       models(:,11) = 3;    % Probability matching
       models(:,13) = 2;    % Lapse
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case {31} % Bisensory standard models with deterministic decision making and lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_dbim';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {41} % Bisensory standard models with deterministic decision making and lapse and simple causal inference (not necessary)
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_dsbim';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) == 4,15) = 3;     % Fixed criterion       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
    case 101 % Bimodal standard models with constant noise
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_bim_const';
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

    case 111 % Bimodal standard models with constant noise without beta
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_2bim_const';
       models(:,11) = 3;    % Probability matching
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case 121 % Bisensory standard models with constant noise without beta and with lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_lbim_const';
       models(:,11) = 3;    % Probability matching
       models(:,13) = 2;    % Lapse
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case 131 % Bisensory standard models with constant noise, deterministic decision making and lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_lbim_const';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case {141} % Bisensory standard models with constant noise, deterministic decision making and lapse and simple causal inference
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_dsbim_const';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) == 4,15) = 3;     % Fixed criterion       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {201} % Bisensory Bayesian models with *constant* noise, deterministic decision making and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Response probability matching
       models = [models; models_pm];       
       options.jobname = 'vest_bayes_cnst';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {202} % Bisensory non-Bayesian models with *constant* noise, deterministic decision making and lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_cnst';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case {203} % Bisensory Bayesian models with *constant* noise, posterior model matching and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayespm_cnst';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,15) = 6; % Model probability matching

   case {204} % Bisensory non-Bayesian models with *constant* noise, deterministic decision making and lapse, NO SIGMA
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_cnst';
       models(:,8) = 3;     % Fixed prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {211} % Bisensory Bayesian models with *cosine* noise, deterministic decision making and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Response probability matching
       models = [models; models_pm];       
       options.jobname = 'vest_bayes';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {212} % Bisensory non-Bayesian models with *cosine* noise, deterministic decision making and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {213} % Bisensory Bayesian models with *cosine* noise, deterministic decision making and lapse (model probability matching)
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayespm_cnst';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,15) = 6; % Model probability matching
       
   case {214} % Bisensory non-Bayesian models with *cosine* noise, deterministic decision making and lapse, NO SIGMA
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_cnst';
       models(:,8) = 3;     % Fixed prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
              
   case {222} % Bisensory non-Bayesian models with *abs-sinusoidal* noise, deterministic decision making and lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = models(:,[1 2]) - 1; % Abs-sinusoidal noise       
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_abs_fix';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   case {232} % Bisensory non-Bayesian models with *cosine* noise, *different w_vis*, deterministic decision making and lapse
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,1) = 3; % Distinct w_vis       
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_wvis_fix';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 6,:) = [];     % Remove Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       
   % CORRELATED PRIOR
       
   case {301} % Bisensory Bayesian models with *constant* noise, deterministic decision making and lapse -- CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_cnst_corr';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_pm = models;
       models_pm(:,15) = 6; % Model probability matching
       models = [models; models_pm];
       models(:,9) = 2; % Fixed delta prior
       models_extra = models;
       models_extra(:,9) = 3;  % Free delta prior
       models = [models; models_extra];
       
   case {304} % Bisensory non-Bayesian models with *constant* noise, deterministic decision making and lapse, NO SIGMA -- CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,9) = 2; % Fixed delta prior
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_cnst_corr';
       models(:,8) = 3;     % Fixed prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 5 | models(:,15) == 6,:) = [];     % Remove Bayesian models and forced fusion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 3;  % Free delta prior
       models = [models; models_extra];

   case {311} % Bisensory Bayesian models with *sinusoidal* noise, deterministic decision making and lapse -- CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_corr';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_pm = models;
       models_pm(:,15) = 6; % Model probability matching
       models = [models; models_pm];
       models(:,9) = 2; % Fixed delta prior
       models_extra = models;
       models_extra(:,9) = 3;  % Free delta prior
       models = [models; models_extra];
       options = setoptions(options,'nstarts',1,1);
       options = setoptions(options,'nsobol',1,1);
       % options.loadinitfromconst = 1;
       options.loadinitfromdisc = 1;
       
   case {314} % Bisensory non-Bayesian models with *sinusoidal* noise, deterministic decision making and lapse, NO SIGMA -- CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,9) = 2; % Fixed delta prior
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_corr';
       models(:,8) = 3;     % Fixed prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 5 | models(:,15) == 6,:) = [];     % Remove Bayesian models and forced fusion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 3;  % Free delta prior
       models = [models; models_extra];
       options = setoptions(options,'nstarts',1,1);
       options = setoptions(options,'nsobol',1,1);
       % options.loadinitfromconst = 1;
       options.loadinitfromdisc = 1;
       
      
    % CORRELATED PRIOR AND DISCRETE STIMULI
    
    case {401} % Bisensory Bayesian models with *constant* noise, deterministic decision making and lapse -- DISCRETE CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_cnst_disc';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_pm = models;
       models_pm(:,15) = 6; % Model probability matching
       models = [models; models_pm];
       models(:,9) = 4; % Fixed delta prior
       models_extra = models;
       models_extra(:,9) = 5;  % Free delta prior
       models = [models; models_extra];
       
   case {403} % Bisensory forced-fusion models with *constant* noise, deterministic decision making and lapse, NO SIGMA -- DISCRETE CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,9) = 4; % Fixed delta prior
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fusion_disc';
       models(:,8) = 5;     % Fixed uniform prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 3 | models(:,15) == 4 | models(:,15) == 6,:) = [];     % Remove Bayesian models and fixed criterion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 5;  % Free delta prior
       models = [models; models_extra];
       
   case {404} % Bisensory non-Bayesian models with *constant* noise, deterministic decision making and lapse, NO SIGMA -- DISCRETE CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,9) = 4; % Fixed delta prior
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_cnst_disc';
       models(:,8) = 5;     % Fixed uniform prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 5 | models(:,15) == 6,:) = [];     % Remove Bayesian models and forced fusion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 5;  % Free delta prior
       models = [models; models_extra];    
    
   case {405} % Bisensory Bayesian models with *constant* noise, deterministic decision making and lapse -- CORRELATED DISCRETE PRIORS ONLY
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models = models(1,:);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_const_disconly';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(:,15) = 7;    % Correlated prior only
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,9) = 4; % Uniform discrete prior on DELTA
       models_extra = models;
       models_extra(:,9) = 5;  % Free discrete Gaussian prior on DELTA
       models = [models; models_extra];
       
   case {406} % Bisensory Bayesian models with *constant* noise, deterministic decision making and lapse -- CORRELATED DISCRETE PRIORS ONLY
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models = models(1,:);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_const_trueprior';
       models(:,8) = 5;     % Uniform discrete prior on SBAR
       models(:,13) = 2;    % Lapse
       models(:,15) = 7;    % Correlated prior only
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,9) = 4; % Uniform discrete prior on DELTA
    
   case {411} % Bisensory Bayesian models with *sinusoidal* noise, deterministic decision making and lapse -- CORRELATED DISCRETE PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_disc';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 2,15) = 1;     % Standard Bayesian
       models(models(:,15) ~= 1 & models(:,15) ~= 2 & models(:,15) ~= 6,:) = [];     % Remove non-Bayesian models       
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_pm = models;
       models_pm(:,15) = 6; % Model probability matching
       models = [models; models_pm];
       models(:,9) = 4; % Uniform discrete prior on DELTA
       models_extra = models;
       models_extra(:,9) = 5;  % Free discrete Gaussian prior on DELTA
       models = [models; models_extra];

   case {413} % Bisensory forced-fusion models with *sinusoidal* noise, deterministic decision making and lapse, NO SIGMA -- DISCRETE CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,9) = 4; % Fixed delta prior
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fusion_disc';
       models(:,8) = 5;     % Fixed uniform prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 3 | models(:,15) == 4 | models(:,15) == 6,:) = [];     % Remove Bayesian models and fixed criterion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 5;  % Free delta prior
       models = [models; models_extra];
       
   case {414} % Bisensory non-Bayesian models with *sinusoidal* noise, deterministic decision making and lapse, NO SIGMA -- DISCRETE CORRELATED PRIORS
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,9) = 4; % Fixed delta prior
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_fix_disc';
       models(:,8) = 5;     % Fixed uniform prior - THESE MODELS DO NOT DEPEND ON SIGMA
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 4,15) = 3;     % Fixed criterion
       models(models(:,15) == 1 | models(:,15) == 2 | models(:,15) == 5 | models(:,15) == 6,:) = [];     % Remove Bayesian models and forced fusion 
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models_extra = models;
       models_extra(:,9) = 5;  % Free delta prior
       models = [models; models_extra];

   case {415} % Bisensory Bayesian models with *sinusoidal* noise, deterministic decision making and lapse -- CORRELATED DISCRETE PRIORS ONLY
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models = models(1,:);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_disconly';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(:,15) = 7;    % Correlated prior only
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,9) = 4; % Uniform discrete prior on DELTA
       models_extra = models;
       models_extra(:,9) = 5;  % Free discrete Gaussian prior on DELTA
       models = [models; models_extra];

   case {416} % Bisensory Bayesian models with *sinusoidal* noise, deterministic decision making and lapse -- CORRELATED DISCRETE PRIORS ONLY
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,models,groupcnd] = VestBMS(options,2,0);
       models = models(1,:);
       models(:,11) = 1;    % BDT
       options.jobname = 'vest_bayes_trueprior';
       models(:,8) = 5;     % Fixed uniform prior on SBAR
       models(:,13) = 2;    % Lapse
       models(:,15) = 7;    % Correlated prior only
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials
       models(:,9) = 4; % Uniform discrete prior on DELTA
       
       
   case {421} % Bisensory probabilistic-fusion models with deterministic decision making and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,~,groupcnd] = VestBMS(options,2,0);
       models = [1 1 1 1, 1 1 1 5, 5 1 1 1, 2 1 9 1; ...
                 1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 9 1; ...
                 5 3 1 1, 1 1 1 5, 5 1 1 1, 2 1 9 1; ...
                 5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 9 1;];
       options.jobname = 'vest_probfusion';
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials       
       
    % LARGE-DISPARITY TRIALS ONLY   
       
    case 501 % Bisensory models

       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_bim_largedisparity';
       dataids(:,2) = setflag(dataids(:,2), [4,5]);     % No categorical trials, no small disparity
        
    case 511 % Bisensory models with constant noise
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       options.jobname = 'vest_bim_const_largedisparity';
       dataids(:,2) = setflag(dataids(:,2), [4,5]);     % No categorical trials, no small disparity
       
%--------------------------------------------------------------------------        
% BISENSORY UNITY JUDGEMENT DATA FITS
       
   case 1001 % Bisensory standard models
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_unity';
       models(:,11) = 3;    % Probability matching
       models(models(:,15) == 5,:) = [];    % Remove forced fusion
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1021 % Bisensory standard models with cosine/constant noise with lapse
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models_const = models;
       models_const(:,[1 2]) = 1; % Constant noise
       models = [models; models_const];
       options.jobname = 'vest_lunity';
       models(:,11) = 3;    % Probability matching
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 5,:) = [];    % Remove forced fusion
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
   case 1031 % Bisensory standard models with random unity judgement
        
       [options,~,groupcnd] = VestBMS(options,2,0);
       options.jobname = 'vest_randunity';
       models = [10 10 1 1, 1 1 1 3, 1 1 1 1, 1 1 5 4 0];
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1041 % Bisensory standard models with cosine/constant noise, BDT and lapse (THIS MIGHT BE WRONG, IT DOES PROBABILITY MATCHING)
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models_const = models;
       models_const(:,[1 2]) = 1; % Constant noise
       models = [models; models_const];
       options.jobname = 'vest_dunity';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 5,:) = [];    % Remove forced fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1051 % Bisensory standard models with cosine/constant noise, probability matching and lapse
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models_const = models;
       models_const(:,[1 2]) = 1; % Constant noise
       models = [models; models_const];
       options.jobname = 'vest_dunity';
       models(:,11) = 3;    % Probability matching
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 5,:) = [];    % Remove forced fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1061 % Bisensory standard models with cosine/constant noise, real BDT and lapse
       [options,models,groupcnd] = VestBMS(options,2,0);
       models_const = models;
       models_const(:,[1 2]) = 1; % Constant noise
       models = [models; models_const];
       options.jobname = 'vest_dbunity';
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
    case 1101 % Bimodal standard models with constant noise
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 3;    % Probability matching
       models(models(:,15) == 5,:) = [];    % Remove forced fusion
       options.jobname = 'vest_unity_const';
       dataids = [(1:11)', zeros(11,1)];       
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
   case 1201 % Bisensory Bayesian models with *constant* noise and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];       
       options.jobname = 'vest_bayes_cnst_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1202 % Bisensory fixed-criterion models with cosine/constant noise, real BDT and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models_const = models;
       models_const(:,[1 2]) = 1; % Constant noise
       models = [models; models_const];
       options.jobname = 'vest_fixed_unity';
       models(:,8) = 3;     % No prior (irrelevant for fixed model)
       models(:,11) = 1;    % BDT
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 1 | models(:,15) == 2,:) = [];    % Remove Bayesian
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
   case 1203 % Bisensory Bayesian models with *cosine* noise and lapse
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];
       options.jobname = 'vest_bayes_ecc_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
   case 1301 % Bisensory Bayesian models with *constant* noise and lapse and CORRELATED priors
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];
       models(:,9) = 3;     % Free correlated prior
       models_prior = models;
       models_prior(:,9) = 2;   % Correlated prior from experiment
       models = [models; models_prior];
       options.jobname = 'vest_bayes_corr_cnst_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials

   case 1311 % Bisensory Bayesian models with *sinusoidal* noise and lapse and CORRELATED priors
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];
       models(:,9) = 3;     % Free correlated prior
       models_prior = models;
       models_prior(:,9) = 2;   % Correlated prior from experiment
       models = [models; models_prior];
       options.jobname = 'vest_bayes_corr_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       options = setoptions(options,'nstarts',1,1);
       options = setoptions(options,'nsobol',1,1);
       %options = setslowoptions(options); % Slow computation       
       % options.loadinitfromconst = 1;
       options.loadinitfromdisc = 1;
       

   case 1401 % Bisensory Bayesian models with *constant* noise and lapse and DISCRETE CORRELATED priors
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,[1 2]) = 1; % Constant noise
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];
       models(:,9) = 5;     % Free correlated prior
       models_prior = models;
       models_prior(:,9) = 4;   % Correlated prior from experiment
       models = [models; models_prior];
       options.jobname = 'vest_bayes_disc_cnst_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
       
   case 1411 % Bisensory Bayesian models with *sinusoidal* noise and lapse and DISCRETE CORRELATED priors
        
       [options,models,groupcnd] = VestBMS(options,2,0);
       models(:,11) = 1;    % BDT
       models_pm = models;
       models_pm(:,11) = 3; % Probability matching
       models = [models; models_pm];
       models(:,9) = 5;     % Free correlated prior
       models_prior = models;
       models_prior(:,9) = 4;   % Correlated prior from experiment
       models = [models; models_prior];
       options.jobname = 'vest_bayes_disc_unity';
       models(:,8) = 2;     % Fixed-mean prior
       models(:,13) = 2;    % Lapse
       models(models(:,15) == 3 | models(:,15) == 4,:) = [];    % Remove fixed criterion
       models(models(:,15) == 5,:) = [];    % Remove probabilistic fusion
       models(:,15) = models(:,15) - 1;     % Remove softness
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
              
%--------------------------------------------------------------------------        
% FULL JOINT DATA FITS

   case 2001 % Full joint standard models (humans)
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_human';       
       models = [ ...
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 2 2 0; ... % Generalized Bayesian
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 4 2 0; ... % Soft fixed criterion
       ];       
       % models(:,11) = 3;    % Probability matching (might require change)
       dataids = [(1:11)', zeros(11,1)];       

   case 2002 % Full joint standard models (monkeys)
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_monkey';
       models = [ ...
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 2 2 0; ... % Generalized Bayesian
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 4 2 0; ... % Soft fixed criterion
       ];       
       % models(:,11) = 3;    % Probability matching (might require change)
       dataids = [12 8; 13 8; 14 8];
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

   case 2011 % Full joint standard models, separate criteria and softmax (humans)
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint2_human';       
       models = [ ...
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 2 3 0; ... % Generalized Bayesian
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 4 3 0; ... % Soft fixed criterion
       ];       
       dataids = [(1:11)', zeros(11,1)];       
       
   case 2021 % Full joint standard models with lapse (humans)
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_human';       
       models = [ ...
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 2 2 0; ... % Generalized Bayesian
           5 3 1 1, 1 1 1 1, 1 1 2 1, 1 1 4 2 0; ... % Soft fixed criterion
       ];       
       models(:,11) = 3;    % Probability matching
       models(:,13) = 2;    % Lapse
       dataids = [(1:11)', zeros(11,1)];

   case 2101 % Full joint standard models with lapse (humans)
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint';
       models = [ ...
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 6 2 0; ... % Bayesian model-PM
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 2 0; ... % Fixed criterion
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
    case 2201   % Full joint standard models with UNCORRELATED priors and ECCENTRIC noise
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_unc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 1 1 0; ... % Bayesian model (same pararameters)
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 1 0; ... % Fixed criterion (same pararameters)
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 1 6 0; ... % Bayesian model (forced fusion)
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 6 0; ... % Fixed criterion (forced fusion)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
    case 2202   % Full joint standard models with UNCORRELATED priors and CONSTANT noise
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_unc_const';
       models = [ ...
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 1 1 0; ... % Bayesian model (same pararameters)
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 1 0; ... % Fixed criterion (same pararameters)
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 1 6 0; ... % Bayesian model (forced fusion)
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 6 0; ... % Fixed criterion (forced fusion)
       ];       
       dataids = [(1:11)', zeros(11,1)];

       
    case 2211   % Full joint standard models with UNCORRELATED/DISCRETE priors and CONSTANT/ECCENTRIC noise and Bayesian probability matching
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_probmatch';
       models = [ ...
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 6 1 0; ... % Eccentric noise, uncorrelated prior (same pararameters)
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 6 1 0; ... % Constant noise, uncorrelated prior (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 6 1 0; ... % Eccentric noise, discrete prior (same pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 6 1 0; ... % Constant noise, discrete prior(same pararameters)
       ];       
       dataids = [(1:11)', zeros(11,1)];

    case 2212   % Full joint standard models with UNCORRELATED/DISCRETE priors and CONSTANT/ECCENTRIC noise and Bayesian model selection
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_modelselect';
       models = [ ...
           5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 8 1 0; ... % Eccentric noise, uncorrelated prior (same pararameters)
           1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 8 1 0; ... % Constant noise, uncorrelated prior (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 8 1 0; ... % Eccentric noise, discrete prior (same pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 8 1 0; ... % Constant noise, discrete prior(same pararameters)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
       
       
    case 2401 % Full joint standard models with discrete priors
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_disc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 1 0; ... % Bayesian model (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 1 0; ... % Fixed criterion (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 5 0; ... % Bayesian model (distinct pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 2 0; ... % Fixed criterion (distinct pararameters)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
   case 2402 % Full joint standard models with discrete priors and forced fusion
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_fusion_disc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 6 0; ... % Bayesian model (forced fusion)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 6 0; ... % Fixed criterion (forced fusion)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
   case 2403 % Full joint standard models with discrete priors and CONSTANT noise
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_disc_const';
       models = [ ...
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 1 0; ... % Bayesian model (same pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 1 0; ... % Fixed criterion (same pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 5 0; ... % Bayesian model (distinct pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 2 0; ... % Fixed criterion (distinct pararameters)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 6 0; ... % Bayesian model (forced fusion)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 6 0; ... % Fixed criterion (forced fusion)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
   case 2404 % Full joint standard models with discrete priors and Bayesian probability matching
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_joint_probmatch_disc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 6 1 0; ... % Bayesian probability matching (eccentric)
           1 1 1 1, 1 1 1 2, 5 1 1 1, 2 1 6 1 0; ... % Bayesian probability matching (constant)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       
       
%--------------------------------------------------------------------------        
% SEMI-JOINT DATA FITS (NO UNISENSORY DATA)
       
   case 3401 % Semi-joint standard models with discrete priors
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 5:7;
       options.jobname = 'vest_semi_disc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 1 0; ... % Bayesian model (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 1 0; ... % Fixed criterion (same pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 5 0; ... % Bayesian model (distinct pararameters)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 2 0; ... % Fixed criterion (distinct pararameters)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 7);     % No unisensory trials
       
   case 3402 % Semi-joint standard models with discrete priors and forced fusion
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
       options = VestBMS(options,2,0);
       groupcnd = 5:7;
       options.jobname = 'vest_semi_fusion_disc';
       models = [ ...
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 1 6 0; ... % Bayesian model (forced fusion)
           5 3 1 1, 1 1 1 2, 5 1 1 1, 2 1 3 6 0; ... % Fixed criterion (forced fusion)
       ];       
       dataids = [(1:11)', zeros(11,1)];
       dataids(:,2) = setflag(dataids(:,2), 7);     % No unisensory trials
       
%--------------------------------------------------------------------------
% UNIMODAL ESTIMATION DATA FITS
% All models have by default: no rescaling, no motor noise, a fixed prior
% (it is irrelevant), no lapse, and softmax decision making.
                
    case {10001}; % Standard unimodal with sinusoidal noise
        
        [options,models,groupcnd] = VestBMS(options,1);
        options.jobname = 'vest_uni';

    case {10101}; % Standard unimodal with constant noise
        
        [options,models,groupcnd] = VestBMS(options,1);
        models(:,[1 2]) = 1; % Constant noise
        options.jobname = 'vest_uni_const';
        
    case 10201 % Standard unimodal with sinusoidal/constant noise and lapse
        
        [options,models,groupcnd] = VestBMS(options,1);
        models_const = models;
        models_const(:,[1 2]) = 1; % Constant noise
        models = [models; models_const];
        models(:,11) = 1;           % BDT (no softmax)
        models(:,13) = 2;           % Lapse        
        options.jobname = 'vest_lapse_uni';        

        
    case 10301 % Standard unimodal with cosine/constant noise and lapse and fixed prior
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
        
        [options,models,groupcnd] = VestBMS(options,1);
        models_const = models;
        models_const(:,[1 2]) = 1; % Constant noise
        models = [models; models_const];
        models(:,8) = 3;            % Fixed prior
        models(:,11) = 1;           % BDT (no softmax)
        models(:,13) = 2;           % Lapse        
        options.jobname = 'vest_noprior_lapse_uni';
        
%--------------------------------------------------------------------------
% MODEL RECOVERY
   
   case 20001 % Bisensory models for the localization task
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       [options,~,groupcnd] = VestBMS(options,2,0);
       models = [5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 6 1; ...    % BPM
            5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 3 1; ...         % CX
            5 3 1 1, 1 1 1 2, 1 1 1 1, 2 1 5 1; ...         % FF
            1 1 1 1, 1 1 1 2, 1 1 1 1, 2 1 6 1];            % BMP-C
       options.jobname = 'vest_modrec_bimloc';
       
        N = 10;          % Ten fake datasets per subject
        Nsubjs = 11;     % Eleven real subjects in total
        dataids = [];    % Take only first five fake datasets per subject
        for i = 1:size(models,1)
            temp = bsxfun(@plus, (1:5)', ((1:Nsubjs)-1)*N) + (i-1)*(N*Nsubjs);
            dataids = [dataids; temp(:)];
        end
        
       dataids = [dataids, zeros(numel(dataids),1)];
       dataids(:,2) = setflag(dataids(:,2), 4);     % No categorical trials

    case 21001 % Bisensory models with unity judgment
        [options,~,groupcnd] = VestBMS(options,2,0);
        models = [5 3 1 1 1 1 1 2 1 1 1 1 2 1 1 1; ... % BP
            5 3 1 1 1 1 1 3 1 1 1 1 2 1 3 1; ...       % CX
            10 10 1 1 1 1 1 3 1 1 1 1 1 1 5 4; ...     % FF
            1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1];          % BP-C       
        options.jobname = 'vest_modrec_unity';

        N = 10;          % Ten fake datasets per subject
        Nsubjs = 11;     % Eleven real subjects in total
        dataids = [];    % Take only first five fake datasets per subject
        for i = 1:size(models,1)
            temp = bsxfun(@plus, (1:5)', ((1:Nsubjs)-1)*N) + (i-1)*(N*Nsubjs);
            dataids = [dataids; temp(:)];
        end
        
       dataids = [dataids, zeros(numel(dataids),1)];
       dataids(:,2) = setflag(dataids(:,2), 3);     % No estimation trials
               
       
    case 22001   % Joint fits model recovery
        % THIS MODEL BELONGS TO THE FINAL MODEL SET
       
       options = VestBMS(options,2,0);
       groupcnd = 1:7;
       options.jobname = 'vest_modrec_joint';
       models = [ ...
           5 3 1 1 1 1 1 2 5 1 1 1 2 1 3 1; ...     % CXD
           5 3 1 1 1 1 1 2 5 1 1 1 2 1 1 1; ...     % BPD
           5 3 1 1 1 1 1 2 1 1 1 1 2 1 1 6; ...     % BPFs
           1 1 1 1 1 1 1 2 1 1 1 1 2 1 3 6; ...     % CXF-CS
           5 3 1 1 1 1 1 2 1 1 1 1 2 1 3 1; ...     % CX
           1 1 1 1 1 1 1 2 5 1 1 1 2 1 3 1 ...      % CXD-C
           ];
       
        Nfake = 10;     % Ten fake datasets per subject
        Nsubjs = 11;    % Eleven real subjects in total
        dataids = [];   % Take only first four fake datasets per subject
        for k = 1:size(models,1)
            temp = bsxfun(@plus, (1:4)', ((1:Nsubjs)-1)*Nfake) + (k-1)*(Nfake*Nsubjs);
            dataids = [dataids; temp(:)];
        end
       dataids = [dataids, zeros(size(dataids,1),1)];
       
end

% Set speed test values
%if any(groupcnd >= 5)
%    options = setoptions(options,'speedtest',10,1);  % Bimodal data
%else
%    options = setoptions(options,'speedtest',0.15,1);  % Unimodal data
%end

% Be verbose
options = setoptions(options,'display','all',1);

% Set conditions
for i = 1:nDatasets; cnd{i} = {groupcnd}; end

% Optimization run
if type > 0 && options.nsamples == 0
    options = setoptions(options,'optfevals',NITER_OPTIMIZATION,1);
    options = setoptions(options,'nsamples',0,1);
    options = setoptions(options,'samplingtemperature',1,1);
    options = setoptions(options,'hessianflag',0,1); % Too expensive - do not compute Hessian
end

return;

    %VESTBMS Define options for Kalpana's experiment
    function [options,models,groupcnd] = VestBMS(options,nStimuli,typeoffset)
        if nargin < 3; typeoffset = []; end
        
        options = setoptions(options,'experimentName','kalpana');
        options = setoptions(options,'bincenters',[-45,-40,-35,-30:2.5:-2.5,-1.25:0.625:1.25,2.5:2.5:30,35,40,45]);
        options = setoptions(options,'respbincenters',[-1,1]);
        options = setoptions(options,'nsamples',NSAMPLES(nStimuli),1);
                
        % Compute bin weights for beliefs about eccentricity-independent likelihoods (unused)
        w = +(abs(options.bincenters) <= 25 & options.bincenters ~= 0);
        w(abs(options.bincenters) < 2.5 & options.bincenters ~= 0) = 0.5;
        for icnd = 1:4
            binweights(:,icnd) = w/sum(w);
        end
        s = -25:5:25; % Mean heading direction
        d = [-40,-20,-10,-5,0,5,10,20,40]'; % Disparity
        s_all = bsxfun(@plus, s, d/2);
        w = zeros(length(options.bincenters),1);
        for iW = 1:length(options.bincenters)
            w(iW) = sum(s_all(:) == options.bincenters(iW));
        end
        for icnd = 5:7
            binweights(:,icnd) = w/sum(w);            
        end
        options = setoptions(options,'binweights',binweights);
        
        % Compute bimodal bin centers
        bincenters_bim{1} = zeros(length(s),length(d));
        bincenters_bim{2} = zeros(length(s),length(d));
        for i1 = 1:length(s)
            for i2 = 1:length(d)
                bincenters_bim{1}(i1,i2) = s(i1) + 0.5*d(i2);
                bincenters_bim{2}(i1,i2) = s(i1) - 0.5*d(i2);
            end
        end
        bincenters_bim{1} = bincenters_bim{1}(:);
        bincenters_bim{2} = bincenters_bim{2}(:);
        bincenters_bim{3} = options.respbincenters;
        options = setoptions(options,'bincenters_bim',bincenters_bim);
        
        if nStimuli == 1 % Unimodal data
            models = standardunimodalmodels;
            groupcnd = 1:4; % All unimodal conditions
            
        else % Bimodal data
            models = standardbimodalmodels;
            modid = type-typeoffset;
            groupcnd = 5:7;
            options = setoptions(options,'unifullthetanumber',modid);
        end
    end

end

function options = setslowoptions(options)
%SETSLOWOPTIONS Set optimization options for slow computations.


end
