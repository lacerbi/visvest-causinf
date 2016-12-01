% CUEBMS_SETUPMODEL initialize and update model parameter structure.
%
% MP = CUEBMS_SETUPMODEL([], THETA, MODEL, INFOSTRUCT) returns a
% model-parameter structure MP initialized with parameter vector THETA 
% according to the model specified in model vector MODEL. INFOSTRUCT is a
% structure whose fields contain additional initialization information.
% THETA can be empty, in which case the parameter structure is only 
% partially initialized.
%
% MP = CUEBMS_SETUPMODEL(MP, THETA) updates existing model-parameter 
% structure MP with global parameter values THETA.
%
% [MP, OUTFLAG] = CUEBMS_SETUPMODEL(...) returns the flag OUTFLAG 
% which takes the value 1 if there were errors, 0 otherwise.

% The MODEL array contains, in order:
%--------------------------------------------------------------------------
% MODEL(1) Visual noise variance:
% 1 Constant (1 param), 2 Quadratic (2 params), 3 Sinusoidal (2 params), 
% 4 Quadratic with single w (2 params), 5 Sinusoidal with single w (2 params)
% 6 Quadratic with single w, base fixed (1 param), 7 Square-abs with single w, base fixed (1 param)
% 8 Quadratic with single-proportional w (2 params), 9 Sinusoidal with single-proportional w (2 params),
% 10 none
%--------------------------------------------------------------------------
% MODEL(2) Vestibular noise variance:
% 1 Constant (1 param), 2 Quadratic (2 params), 3 Sinusoidal (2 params), 
% 4 Quadratic with w equal to w_vis (1 param), 5 Sinusoidal with w equal to w_vis (1 param).
% 6 Quadratic with base fixed (1 param), 7 Sinusoidal with base fixed (1 param)
% 8 Quadratic with base fixed and w equal to w_vis, 9 Sinusoidal with base fixed and w equal to w_vis
% 10 none
%--------------------------------------------------------------------------
% MODEL(3) Sensory noise source:
% 1 None, 2 From unimodal data, 3 Proportional from data (1 param), 
% 4 Vestibular base proportional from data (1 param), 5 Vestibular-video base independently proportional from data (2 params), 
% 6 Vestibular fully proportional from data (1 param), 7 Vestibular-video independently fully proportional from data (2 params), 
% 8 Vestibular-video independently fully proportional from data, multiple video (n+1 params),
% 9 Vestibular-video base independently proportional from data, leave w parameter(s) free (2 params)
%--------------------------------------------------------------------------
% MODEL(4) Visual noise subjective variance (unused):
% 1 Equal to external.
%--------------------------------------------------------------------------
% MODEL(5) Vestibular noise subjective variance (unused):
% 1 Equal to external.
%--------------------------------------------------------------------------
% MODEL(6) Visual likelihood additional parameters (unused):
% 1 None, 2 Rescaling (1 param), 3 Rescaling and shift (2 params), 4 From unimodal data.
%--------------------------------------------------------------------------
% MODEL(7) Vestibular likelihood additional parameters (unused):
% 1 None, 2 Rescaling (1 param), 3 Rescaling and shift (2 params), 4 From unimodal data.
%--------------------------------------------------------------------------
% MODEL(8) Prior model:
% 1 Free Gaussian (2 params), 2 Centered Gaussian (1 param), 3 Fixed Gaussian
% 4 From unimodal data, 5 Fixed uniform
%--------------------------------------------------------------------------
% MODEL(9) Prior model (correlated or not):
% 1 Uncorrelated prior, 
%  2 Correlated prior (fixed), 3 Correlated prior (free), 
%  4 Discrete correlated prior (fixed), 5 Discrete correlated prior (free)
%--------------------------------------------------------------------------
% MODEL(10) Loss model (unused):
% 1 Unused
%--------------------------------------------------------------------------
% MODEL(11) Decision making model:
% 1 BDT, 2 Power-law (1 param), 3 Probability matching
%--------------------------------------------------------------------------
% MODEL(12) Unused
% 1 Unused
%--------------------------------------------------------------------------
% MODEL(13) Lapse model:
% 1 No lapse, 2 Lapse (1 param), 3 From data.
%--------------------------------------------------------------------------
% MODEL(14) Causal inference model (Bimodal only, unused):
% 1 Unused
%--------------------------------------------------------------------------
% MODEL(15) Model weighting model (Bimodal only):
% 1 Bayesian posterior (1 param), Generalized Bayesian posterior (2 params),
% 3 Criterion on x (1 param), 4 Soft criterion on x (2 params),
% 5 Forced fusion, 6 Bayesian posterior probability matching (1 param)
% 7 Bayesian posterior, true p_common, 8 Bayesian posterior model selection (1 param)
% 9 Probabilistic fusion
%--------------------------------------------------------------------------
% MODEL(16) Report of unity model (Bimodal only):
% 1 Standard, 2 Separate criterion parameter (1 param), 
% 3 Separate criterion and gamma (2 params), 4 Free probability (3 params),
% 5 Separate criterion and sigmadelta (2 params),
% 6 Forced fusion on localization

function [mp,exitflag] = VestBMS_setupModel(mp,theta,model,infostruct)

if ~exist('theta', 'var'); theta = []; end
if ~exist('model', 'var'); model = []; end
if ~exist('infostruct', 'var'); infostruct = []; end

% Initialize model
if isempty(mp)
    [mp, exitflag] = initModel(model, infostruct);
    if ~isempty(theta) && exitflag == 0
        [mp, exitflag] = updateModel(mp,theta);
    end
else
    [mp, exitflag] = updateModel(mp,theta);
end

end

% INITMODEL Initialize model.
function [mp, outflag] = initModel(model, infostruct)
    mp = [];
    outflag = 0;
    
    % Take information from infostruct
    cnd = infostruct.cnd;
                    
    % Number of function calls in the optimization
    mp.funcalls = 0;
    
    mp.model = model;
    mp.cnd = cnd;
    mp.ncnd = length(cnd);
        
    % For advanced sampling parameters - is kappa going to start from zero?
    mp.kappaminusoneflag = 1;
    
    % Multiplier to slice sampling step size
    if ~isfield(infostruct, 'samplingtemperature'); infostruct.samplingtemperature = []; end
    if isempty(infostruct.samplingtemperature); infostruct.samplingtemperature = 1; end
    mp.samplingtemperature = infostruct.samplingtemperature;
       
    % Parameters from unimodal conditions
    if ~isfield(infostruct, 'unifulltheta'); infostruct.unifulltheta = []; end
    mp.unifulltheta = infostruct.unifulltheta;
    
    % Binned intervals
    mp.bincenters = infostruct.bincenters;
    
    % Bin weights (for computing eccentricity-independent likelihood; unused)
    if ~isfield(infostruct, 'binweights'); infostruct.binweights = []; end
    mp.binweights = infostruct.binweights;
    if ~isempty(mp.bincenters) && isempty(mp.binweights)
        mp.binweights = ones(length(mp.bincenters),mp.ncnd)/length(mp.bincenters);
    end
    
    % Maximum stimulus range
    mp.MAXRNG = infostruct.MAXRNG; 
    if isfield(infostruct,'MAXDELTA'); mp.MAXDELTA = infostruct.MAXDELTA; end
    
    % Fixed parameters
    % if ~isfield(infostruct, 'fixedtheta'); infostruct.fixedtheta = []; end
    % mp.fixedtheta = infostruct.fixedtheta;
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup model parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Bounds are maximal lower/upper bounds, reasonable lower/upper bounds,
    % and reasonable width of the posterior (for slice sampling)

    pbounds = [];
    params = [];
    
    % Number of parameters
    mp.nparams = zeros(1, 16);
    for i = 1:size(mp.nparams, 2); pbounds{i} = []; params{i} = []; end
    
    % Parameter shared across conditions
    if ~isfield(infostruct, 'sharedparams'); infostruct.sharedparams = []; end
    if isempty(infostruct.sharedparams);
       % By default, all parameters are shared
       infostruct.sharedparams = ones(1, size(mp.nparams, 2));
    end
    mp.sharedparams = infostruct.sharedparams;
    %if ~mp.sharedparams(5)
    %    error('Sensory mapping parameters need to be shared between sessions.');
    %end
    
    mp.theta = [];
    
    sigmabound = log([0.5 80, 1 40, exp(2)]);
    wbound = [0 1, 0 0.5, 0.15];
    softmax_bounds = log([1e-4 2, 0.01 0.5, exp(1)]);

    % Visual noise variance
    % Set default values
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.sigma_vis_low = NaN;
        mp.fulltheta{icnd}.sigma_vis_med = NaN;
        mp.fulltheta{icnd}.sigma_vis_high = NaN;
        mp.fulltheta{icnd}.w_vis_low = 0;
        mp.fulltheta{icnd}.w_vis_med = 0;
        mp.fulltheta{icnd}.w_vis_high = 0;
    end
    switch model(1)
        case 1 % Constant variance (log scale)
        temp = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high'};
        for iNoise = 1:3            
            if any(mp.cnd == iNoise) || any(mp.cnd == (iNoise+4))
                mp.nparams(1) = mp.nparams(1) + 1;
                params{1}{end+1} = temp{iNoise};
            end
        end
        pbounds{1} = repmat(sigmabound, [mp.nparams(1) 1]);

        case {2, 3} % Quadratic or sinusoidal variance (log scale)
        temp = {'sigma_vis_low'; 'w_vis_low'; 'sigma_vis_med'; ...
                'w_vis_med'; 'sigma_vis_high'; 'w_vis_high'};
        for iNoise = 1:3
            if any(mp.cnd == iNoise) || any(mp.cnd == (iNoise+4))
                mp.nparams(1) = mp.nparams(1) + 2;
                params{1}{end+1} = temp{iNoise*2-1};
                params{1}{end+1} = temp{iNoise*2};
            end
        end
        pbounds{1} = repmat([sigmabound; wbound], [mp.nparams(1)/2, 1]);

        case {4,5,8,9} % Quadratic or sinusoidal variance with single visual w (log scale)
        temp = {'sigma_vis_low'; 'sigma_vis_med'; 'sigma_vis_high'; 'w_vis'};
        for iNoise = 1:3
            if any(mp.cnd == iNoise) || any(mp.cnd == (iNoise+4))
                mp.nparams(1) = mp.nparams(1) + 1;
                params{1}{end+1} = temp{iNoise};
            end
        end
        mp.nparams(1) = mp.nparams(1) + 1;
        params{1}{end+1} = temp{end};
        pbounds{1} = [repmat(sigmabound, [mp.nparams(1)-1, 1]); wbound];
        
        case {6,7} % Quadratic or sinusoidal variance with fixed sigma (log scale)
        mp.nparams(1) = 1;
        params{1}{1} = 'w_vis';
        pbounds{1} = wbound;
        
        case 10 % None

        otherwise
            error('Unsupported sensory variance model.');
    end

    % Vestibular noise variance
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.sigma_vest = 0;
        mp.fulltheta{icnd}.w_vest = 0;
    end
    switch model(2)
        case 1 % Constant variance (log scale)
            mp.nparams(2) = 1;
            pbounds{2} = sigmabound;
            params{2} = {'sigma_vest'};
        case {2, 3} % Quadratic or sinusoidal variance (log scale)
            mp.nparams(2) = 2;
            pbounds{2} = [sigmabound; wbound];
            params{2} = {'sigma_vest'; 'w_vest'};                        
        case {4, 5} % Quadratic of sinusoidal variance with w fixed to w_vis (log scale)
            mp.nparams(2) = 1;
            pbounds{2} = sigmabound;
            params{2} = {'sigma_vest'};
        case {6, 7} % Quadratic or sinusoidal variance with base fixed (log scale)
            mp.nparams(2) = 1;
            pbounds{2} = wbound;
            params{2} = {'w_vest'};
        case {8, 9} % Quadratic or sinusoidal variance with base fixed and w equal to w_vis (log scale)
            mp.nparams(2) = 0;
        case 10 % None            
        otherwise; error('Unsupported sensory variance model.');
    end
    
    % Sensory noise source
    if model(1) == 1 && model(2) == 1 || model(3) == 9
        paramstring = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest'};
    elseif (model(1) == 2 || model(1) == 3) && any(model(2) == [2 3])
        paramstring = {'sigma_vis_low', 'w_vis_low', 'sigma_vis_med', ...
                    'w_vis_med', 'sigma_vis_high', 'w_vis_high', 'sigma_vest', 'w_vest'};
    elseif any(model(1) == 4:9) && any(model(2) == [2 3])
        paramstring = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', ...
            'w_vis', 'sigma_vest', 'w_vest'};
    elseif any(model(1) == 4:9) && any(model(2) == [4 5])
        paramstring = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', ...
            'w_vis', 'sigma_vest'};
    elseif model(3) > 1
        error('Unsupported noise combination.')
    end
    switch model(3)
        case 1 % Standard
        case 2 % From unimodal data
            mp.nparams(1:2) = [0 0];
            params{1} = []; params{2} = [];
            pbounds{1} = []; pbounds{2} = [];
            readParametersFromUnimodalData(paramstring);            
        case {3, 4, 6} % Proportional to unimodal data (1 param)
            mp.nparams(1:3) = [0 0 1];
            params{1} = []; params{2} = [];
            pbounds{1} = []; pbounds{2} = [];
            params{3} = {'k_sens'};
            pbounds{3} = [0.1 6, 0.5 1.5, 0.2];
            readParametersFromUnimodalData(paramstring);
        case {5, 7, 9} % Independently proportional to unimodal data (2 params; model(3) = 9 has free w)            
            if model(3) == 9
                mp.nparams(3) = 2; % Keep sensory parameters
            else
                mp.nparams(1:3) = [0 0 2]; % Remove all sensory parameters
                params{1} = []; params{2} = [];
                pbounds{1} = []; pbounds{2} = [];
            end
            params{3} = {'k_vis', 'k_vest'};
            pbounds{3} = [0.1 6, 0.5 1.5, 0.2; 0.1 4, 0.5 1.5, 0.2];
            readParametersFromUnimodalData(paramstring);            
        case 8 % Independently proportional to unimodal data (1+n params)
            np = 1 + mp.ncnd;
            mp.nparams(1:3) = [0 0 np];
            params{1} = []; params{2} = [];
            pbounds{1} = []; pbounds{2} = [];
            pbounds{3} = []; params{3} = [];
            for i = 1:mp.ncnd
                params{3}{end+1} = ['k_vis#' num2str(i)];
                pbounds{3} = [pbounds{3}; 0.1 10, 0.5 1.5, 0.2];
            end
            params{3}{end+1} = 'k_vest';
            pbounds{3} = [pbounds{3}; 0.1 5, 0.5 1.5, 0.2];
            readParametersFromUnimodalData(paramstring);            
        otherwise; error('Unsupported sensory noise source.');
    end
    
    % Visual subjective likelihood variance (unused)
    % Set default values
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.sigmalike_vis_low = NaN;
        mp.fulltheta{icnd}.sigmalike_vis_med = NaN;
        mp.fulltheta{icnd}.sigmalike_vis_high = NaN;
        mp.fulltheta{icnd}.wlike_vis_low = 0;
        mp.fulltheta{icnd}.wlike_vis_med = 0;
        mp.fulltheta{icnd}.wlike_vis_high = 0;
    end
    switch model(4)
        case 1 % Same as external
        otherwise; error('Unsupported internal sensory likelihood model.');
    end

    % Vestibular subjective likelihood variance (unused)
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.sigmalike_vest = 0;
        mp.fulltheta{icnd}.wlike_vest = 0;
    end
    switch model(5)
        case 1 % Same as external
        otherwise; error('Unsupported internal sensory likelihood model.');
    end

    % Visual and vestibular additional noise parameters (unused)
    string = {'_vis', '_vest'};
    for i = 6:7
        switch model(i)
            case 1 % None                
                for icnd = 1:mp.ncnd
                    mp.fulltheta{icnd}.(['alpha' string{i-5}]) = 1;
                    mp.fulltheta{icnd}.(['beta' string{i-5}]) = 0;
                end
            case 2 % Rescaling
                mp.nparams(i) = 1;
                pbounds{i} = [0.25 3, 0.5 1.5, 0.3];
                params{i} = {['alpha' string{i-5}]};                        
                for icnd = 1:mp.ncnd
                    mp.fulltheta{icnd}.(['beta' string{i-5}]) = 0;
                end                
            case 3 % Rescaling and shift
                mp.nparams(i) = 2;
                pbounds{i} = [0.25 3, 0.5 1.5, 0.3; -20 20, -3 3, 4];
                params{i} = {['alpha' string{i-5}], ['beta' string{i-5}]};                        
            case 4 % From unimodal data
                readParametersFromUnimodalData({['alpha' string{i-5}], ['beta' string{i-5}]});
            otherwise; error('Unsupported additional likelihood noise model.');
        end
    end
    
    %----------------------------------------------------------------------    
    % Prior model
    % Add model for separate prior between uni-sensory and bi-sensory if
    % you ever fit those together
    
    % Set default values    
    for icnd = 1:mp.ncnd
        % Default prior in the experiment
        mp.fulltheta{icnd}.priormu = 0; % Empirical prior mean
        mp.fulltheta{icnd}.priorsigma = 15.8114; % Empirical prior SD
    end
    priormu_bounds = [-mp.MAXRNG mp.MAXRNG, -5 5, 5];
    priorsigma_bounds = log([1 2*mp.MAXRNG, 4 mp.MAXRNG, exp(1)]);
    
    switch model(8)
        case 1 % Gaussian (2-params)
            mp.nparams(8) = 2;
            pbounds{8} = [priormu_bounds; priorsigma_bounds];
            params{8} = {'priormu', 'priorsigma'};
        case 2 % Centered Gaussian (1 params)
            mp.nparams(8) = 1;
            pbounds{8} = priorsigma_bounds;
            params{8} = {'priorsigma'};
        case 3 % Fixed Gaussian
        case 4 % From unimodal data
            readParametersFromUnimodalData({'priormu', 'priorsigma'}); 
        case 5 % Fixed uniform
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.priorsigma = Inf; end
        otherwise; error('Unsupported prior model.');
    end
        
    % Prior model, part two
    switch model(9)
        case 1 % Standard (uncorrelated prior)
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.priorsigmadelta = 0; end
        case 2 % Fixed Gaussian prior on Delta (correlated prior) 
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.priorsigmadelta = 21.7307; end
        case 3 % Free Gaussian prior on Delta (correlated prior)
            mp.nparams(9) = 1;
            pbounds{9} = priorsigma_bounds;
            params{9} = {'priorsigmadelta'};
        case 4 % Uniform prior on discrete stimuli (discrete correlated prior) 
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.priorsigmadelta = Inf; end
        case 5 % Free Gaussian prior on discrete stimuli (discrete correlated prior)
            mp.nparams(9) = 1;
            pbounds{9} = priorsigma_bounds;
            params{9} = {'priorsigmadelta'};            
        otherwise; error('Unsupported prior model (correlation).');        
    end
    
    % Loss function model (unused)
    switch model(10)
        case 1 % Unused
        otherwise; error('Unsupported model.');        
    end
        
    % Decision making model
    switch model(11)
        case 1 % Standard BDT
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.tau_softmax = 0; end
        case 2 % Power-function
            mp.nparams(11) = 1;
            pbounds{11} = softmax_bounds;
            params{11} = {'tau_softmax'}; % Temperature parameter
        case 3 % Probability matching
            for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.tau_softmax = 1; end
        otherwise; error('Unsupported decision making model.');
    end
    
    % Motor noise model (unused)
    switch model(12)
        case 1 % Unused
        otherwise; error('Unsupported model.');        
    end
        
    % Lapse model (unused)
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.lambda = 0; 
    end    
    switch model(13)
        case 1 % No lapse
        case 2 % Lapse (1-param)
            mp.nparams(13) = 1;
            pbounds{13} = [0 1, 0.01 0.2, 0.04];
            params{13} = {'lambda'};
        case 3 % From unimodal data
            readParametersFromUnimodalData({'lambda', 'lambdasigma', 'sigmalapse_vis', 'sigmalapse_vest'});      
        otherwise; error('Unsupported lapse model.');
    end
    
    % Only for bimodal datasets
    if any(mp.cnd > 4)
        
        kcommon_bounds = log([0.25 mp.MAXRNG*2, 1 mp.MAXRNG/2, exp(1)]);
        
        % Causal inference model (Bimodal only)
        for icnd = 1:mp.ncnd
            mp.fulltheta{icnd}.pcommon = 0.5;
            mp.fulltheta{icnd}.kcommon = 0;
            mp.fulltheta{icnd}.pcommon_unity = 0.5;
            mp.fulltheta{icnd}.kcommon_unity = 0;
        end
        
        switch model(14)
            case 1 % Unused
            otherwise; error('Unsupported causal inference model.');
        end

        switch model(15)
            case {1,6,8,9} % Bayesian posterior or probabilistic fusion
                mp.nparams(15) = 1;
                pbounds{15} = [0 1, 0.1 0.9, 0.1];
                params{15} = {'pcommon'};
            case 2 % Generalized Bayesian weight
                mp.nparams(15) = 2;
                pbounds{15} = [0 1, 0.1 0.9, 0.1; softmax_bounds];
                params{15} = {'pcommon','invgamma_causinf'};                
            case 3 % Criterion on x (1-param)
                mp.nparams(15) = 1;
                pbounds{15} = kcommon_bounds;
                params{15} = {'kcommon'};                               
            case 4 % Soft criterion on x (2-params)
                mp.nparams(15) = 2;
                pbounds{15} = [kcommon_bounds; kcommon_bounds];
                params{15} = {'kcommon','tau_causinf'};
            case 5 % Forced fusion
                for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.pcommon = 1; end
            case 7 % Bayesian posterior, true pcommon is 1/5 (checked with data)
                for icnd = 1:mp.ncnd; mp.fulltheta{icnd}.pcommon = 0.2; end
        end
        % For these models criteria are not shared across conditions
        %if model(15) == 3 || model(15) == 5 || model(15) == 7
        %    mp.sharedparams(15) = 0;
        %end

        % Report of unity criterion
        switch model(16)
            case 1 % Standard
            case 2 % Separate criterion (based on 15)                
                switch model(15)
                    case {1,2,6,7,8,9} % Bayesian or probabilistic fusion (1-param)
                        mp.nparams(16) = 1;
                        pbounds{16} = [0 1, 0.1 0.9, 0.1];
                        params{16} = {'pcommon_unity'};
                    case {3,4} % Criterion on x (1-param)
                        mp.nparams(16) = 1;
                        pbounds{16} = kcommon_bounds;
                        params{16} = {'kcommon_unity'};
                end
            case 3 % Separate criterion and gamma (based on 15)
                switch model(15)
                    case {1,2,6,7,8,9} % Bayesian (2-params)
                        mp.nparams(16) = 2;
                        pbounds{16} = [0 1, 0.1 0.9, 0.1; softmax_bounds];
                        params{16} = {'pcommon_unity','invgamma_causinf_unity'};
                    case {3,4} % Criterion on x (2-params)
                        mp.nparams(16) = 2;
                        pbounds{16} = [kcommon_bounds; kcommon_bounds];
                        params{16} = {'kcommon_unity','invgamma_causinf_unity'};
                end
            case 4
                mp.nparams(16) = 3;
                pbounds{16} = [0 1, 0.1 0.9, 0.1; 0 1, 0.1 0.9, 0.1; 0 1, 0.1 0.9, 0.1];
                params{16} = {'random_unity_low','random_unity_med','random_unity_high'};
            case 5  % Separate criterion and priorsigmadelta
                switch model(15)
                    case {1,2,6,7,8,9} % Bayesian (2-params)
                        mp.nparams(16) = 2;
                        pbounds{16} = [0 1, 0.1 0.9, 0.1; priorsigma_bounds];
                        params{16} = {'pcommon_unity','priorsigmadelta_unity'};
                    case {3,4} % Criterion on x (2-params)
                        mp.nparams(16) = 2;
                        pbounds{16} = [kcommon_bounds; priorsigma_bounds];
                        params{16} = {'kcommon_unity','priorsigmadelta_unity'};
                end
            case 6 % Forced fusion on localization                
            otherwise; error('Unsupported report of unity criterion model.');
        end
        
    end
    
    %----------------------------------------------------------------------
    % Finalize parameters

    % Unwrap parameters
    pboundslist = []; paramslist = [];
    for i = 1:size(mp.nparams, 2)
        if ~isempty(pbounds{i})
            pboundslist = [pboundslist; pbounds{i}];
            for j = 1:length(params{i}); paramslist{end+1} = params{i}{j}; end
        end
    end

    mp.nparams = mp.nparams + mp.nparams.*(1-mp.sharedparams)*(mp.ncnd-1);
    
    for icnd = 2:mp.ncnd
        for ii = 1:size(mp.sharedparams, 2)
            if ~mp.sharedparams(ii)
                pboundslist = [pboundslist; pbounds{ii}];
                for j = 1:length(params{ii}); paramslist{end+1} = [params{ii}{j} '#' num2str(icnd)]; end
            end
        end
    end
        
    % Total number of parameters
    mp.ntotparams = sum(mp.nparams, 2);
    
    % for icnd = 1:mp.cnd; mp.fulltheta{icnd} = zeros(1, mp.nparams); end
    
    % Parameters extreme lower bounds and upper bounds
    mp.bounds.LB = zeros(1, mp.ntotparams); mp.bounds.UB = zeros(1, mp.ntotparams); 
    
    % Parameters reasonable lower bounds and upper bounds
    mp.bounds.RLB = zeros(1, mp.ntotparams); mp.bounds.RUB = zeros(1, mp.ntotparams);
    
    % Rescale parameter by 'temperature' parameter
    if isscalar(mp.samplingtemperature)
        mp.samplingtemperature = mp.samplingtemperature*ones(1, mp.ntotparams);
    end
    
    for k = 1:mp.ntotparams
        mp.bounds.LB(k) = pboundslist(k, 1);
        mp.bounds.UB(k) = pboundslist(k, 2);
        mp.bounds.RLB(k) = pboundslist(k, 3);
        mp.bounds.RUB(k) = pboundslist(k, 4);
        mp.bounds.SCALE(k) = pboundslist(k, 5)*mp.samplingtemperature(k);                
    end    
    
    mp.paramstree = params;
    mp.params = paramslist;
    
    return;
    
    function readParametersFromUnimodalData(readparams)    
        for iicnd = 1:mp.ncnd
            for iiParam = 1:length(readparams)
                unicnd = mp.cnd(iicnd) - 4; % Look at the same level of noise
                mp.fulltheta{iicnd}.(readparams{iiParam}) = mp.unifulltheta{unicnd}.(readparams{iiParam});
            end
        end                                        
    end
end

% UPDATEMODEL Adjust global model parameters based on current
% model and given parameter vector. Returns a flag OUTFLAG which is
% one if everything went okay or zero otherwise.
function [mp,exitflag] = updateModel(mp,theta)

    mp.theta = theta;
    model = mp.model;

    % Set XGRID values for different types of computation (precision)
    if ~isfield(mp, 'computation') || isempty(mp.computation)
        mp.computation = 1;  % Standard level of precision
    end
    
    if ischar(mp.computation)
        switch mp.computation
            case 'coarse'; mp.computation = 0; % Coarse grid
            case {'normal','hessian'}; mp.computation = 1; % Fine grid
            case 'precise'; mp.computation = 2; % Ultra-fine grid
        end
    end
    
    if model(9) == 1 % Uncorrelated priors, everything is easy
        gridPoints = [  501 501 501 501, 401 401 401; ...       % Coarse
                        501 501 501 501, 401 401 401; ...       % Fine
                        501 501 501 501, 401 401 401];          % Ultra-fine
   
                    % Correlated priors but with constant noise or discrete stimuli
    elseif ( model(9) > 1 && model(1) == 1 && model(2) == 1 ) || ...
            ( model(9) == 4 || model(9) == 5 )
        gridPoints = [  501 501 501 501, 401 401 401; ...       % Coarse
                        501 501 501 501, 401 401 401; ...       % Fine
                        501 501 501 501, 401 401 401];          % Ultra-fine
    else                    % Correlated prior with eccentric noise
        nslow = 301;
        gridPoints = [  501 501 501 501, nslow nslow nslow; ...       % Coarse
                        501 501 501 501, nslow nslow nslow; ...       % Fine
                        501 501 501 501, nslow nslow nslow];          % Ultra-fine
    end
        
    if mp.computation == round(mp.computation)  % Check if integer
        mp.XGRID = gridPoints(mp.computation+1,:);
    else
        mp.XGRID = round(interp1([0 1 2], gridPoints, mp.computation, 'pchip'));
        f = mod(mp.XGRID,2) == 0;
        mp.XGRID(f) = mp.XGRID(f) + 1; % XGRID is always a odd integer number 
    end
    
    exitflag = 0;
    actparam = 1;    

   % Information for log priors 

    for icnd = 1:mp.ncnd
        thiscnd = mp.cnd(icnd);

        % Visual noise variance    
        switch model(1)
            case 1 % Constant variance (log scale)
                switch mp.nparams(1)
                    case 1; updateparams(icnd, 1, {'exp'});
                    case 2; updateparams(icnd, 1, {'exp', 'exp'});
                    case 3; updateparams(icnd, 1, {'exp', 'exp', 'exp'});
                end

            case {2,3} % Quadratic or sinusoidal variance (log scale)
                switch mp.nparams(1)
                    case 2; updateparams(icnd, 1, {'exp', 'id'});
                    case 4; updateparams(icnd, 1, {'exp', 'id', 'exp', 'id'});
                    case 6; updateparams(icnd, 1, {'exp', 'id', 'exp', 'id', 'exp', 'id'});
                end
                
            case {4,5} % Quadratic or sinusoidal variance with single visual w (log scale)                
                switch mp.nparams(1)
                    case 2; updateparams(icnd, 1, {'exp', 'id'});
                    case 3; updateparams(icnd, 1, {'exp', 'exp', 'id'});
                    case 4; updateparams(icnd, 1, {'exp', 'exp', 'exp', 'id'});
                end
                mp.fulltheta{icnd}.w_vis_low = mp.fulltheta{icnd}.w_vis;
                mp.fulltheta{icnd}.w_vis_med = mp.fulltheta{icnd}.w_vis;
                mp.fulltheta{icnd}.w_vis_high = mp.fulltheta{icnd}.w_vis;

            case {6,7} % Quadratic or sinusoidal variance with single visual w and base fixed (log scale)                
                updateparams(icnd, 1, {'id'});
                mp.fulltheta{icnd}.w_vis_low = mp.fulltheta{icnd}.w_vis;
                mp.fulltheta{icnd}.w_vis_med = mp.fulltheta{icnd}.w_vis;
                mp.fulltheta{icnd}.w_vis_high = mp.fulltheta{icnd}.w_vis;

            case {8,9} % Quadratic or sinusoidal variance with single-proportional visual w (log scale)                
                switch mp.nparams(1)
                    case 2; updateparams(icnd, 1, {'exp', 'id'});
                    case 3; updateparams(icnd, 1, {'exp', 'exp', 'id'});
                    case 4; updateparams(icnd, 1, {'exp', 'exp', 'exp', 'id'});
                end
                mp.fulltheta{icnd}.w_vis_low = mp.fulltheta{icnd}.w_vis*mp.fulltheta{icnd}.sigma_vis_low/mp.fulltheta{icnd}.sigma_vis_med;
                mp.fulltheta{icnd}.w_vis_med = mp.fulltheta{icnd}.w_vis;
                mp.fulltheta{icnd}.w_vis_high = mp.fulltheta{icnd}.w_vis*mp.fulltheta{icnd}.sigma_vis_high/mp.fulltheta{icnd}.sigma_vis_med;
        end
        
        % Vestibular noise variance
        duplicate_wvis = 0;
        switch model(2)
            case 1 % Constant variance (log scale)
                updateparams(icnd, 2, {'exp'});
            case {2, 3} % Quadratic or square-abs variance (log scale)
                updateparams(icnd, 2, {'exp', 'id'});
            case {4, 5} % Quadratic or square-abs with fixed w equal to w_vis (log scale)
                updateparams(icnd, 2, {'exp'});
                duplicate_wvis = 1;
            case {6, 7} % Quadratic or square-abs variance with base fixed
                updateparams(icnd, 2, {'id'});
            case {8, 9} % Quadratic or square-abs variance with base fixed and w equal to w_vis
                duplicate_wvis = 1;
        end
        % For these conditions w_vest is equal to (average) w_vis
        if duplicate_wvis
            % Check if the conditions exist
            if isnan(mp.fulltheta{icnd}.sigma_vis_high) && isnan(mp.fulltheta{icnd}.sigma_vis_med)
                mp.fulltheta{icnd}.w_vest = mp.fulltheta{icnd}.w_vis_low;
            elseif isnan(mp.fulltheta{icnd}.sigma_vis_high)
                mp.fulltheta{icnd}.w_vest = sqrt(mean([mp.fulltheta{icnd}.w_vis_low.^2, mp.fulltheta{icnd}.w_vis_med.^2]));                       
            else
                mp.fulltheta{icnd}.w_vest = sqrt(mean([mp.fulltheta{icnd}.w_vis_low.^2, mp.fulltheta{icnd}.w_vis_med.^2, mp.fulltheta{icnd}.w_vis_high.^2]));
            end
        end
        
        % Sensory noise source
        switch model(3)
            case 1 % Standard
            case 2 % Full from unimodal data
            case {3, 4} % Proportional to unimodal data (1 param)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    if model(3) == 3 % All noise is proportional
                        for iParam = 1:size(readparams, 2)
                            mp.fulltheta{icnd}.(readparams{iParam}) = theta(actparam)*mp.unifulltheta{unicnd}.(readparams{iParam});
                        end
                    elseif model(3) == 4 % Only vestibular noise is proportional
                        mp.fulltheta{icnd}.sigma_vest = theta(actparam)*mp.unifulltheta{unicnd}.sigma_vest;
                    end
                    mp.fulltheta{icnd}.k_sens = theta(actparam);
                    actparam = actparam + 1;
                end
            case 5 % Independently proportional to unimodal data (2 params)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    for iParam = 1:3
                        mp.fulltheta{icnd}.(readparams{iParam}) = theta(actparam)*mp.unifulltheta{unicnd}.(readparams{iParam});
                    end
                    mp.fulltheta{icnd}.sigma_vest = theta(actparam+1)*mp.unifulltheta{unicnd}.sigma_vest;
                    mp.fulltheta{icnd}.k_vis = theta(actparam);
                    mp.fulltheta{icnd}.k_vest = theta(actparam+1);
                    actparam = actparam + 2;
                end
            case 6 % Vestibular noise fully proportional to unimodal data (1 param)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest', ...
                    'w_vis_low', 'w_vis_med', 'w_vis_high', 'w_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    mp.fulltheta{icnd}.sigma_vest = theta(actparam)*mp.unifulltheta{unicnd}.sigma_vest;
                    mp.fulltheta{icnd}.w_vest = theta(actparam)*mp.unifulltheta{unicnd}.w_vest;
                    mp.fulltheta{icnd}.k_sens = theta(actparam);
                    actparam = actparam + 1;
                end
            case 7 % Independently fully proportional to unimodal data (2 params)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest', ...
                    'w_vis_low', 'w_vis_med', 'w_vis_high', 'w_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    for iParam = 1:3
                        mp.fulltheta{icnd}.(readparams{iParam}) = theta(actparam)*mp.unifulltheta{unicnd}.(readparams{iParam});
                        mp.fulltheta{icnd}.(readparams{iParam+4}) = theta(actparam)*mp.unifulltheta{unicnd}.(readparams{iParam+4});
                    end
                    mp.fulltheta{icnd}.sigma_vest = theta(actparam+1)*mp.unifulltheta{unicnd}.sigma_vest;
                    mp.fulltheta{icnd}.w_vest = theta(actparam+1)*mp.unifulltheta{unicnd}.w_vest;
                    mp.fulltheta{icnd}.k_vis = theta(actparam);
                    mp.fulltheta{icnd}.k_vest = theta(actparam+1);
                    actparam = actparam + 2;
                end
            case 8 % Independently and multiply fully proportional to unimodal data (1+ncnd params)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest', ...
                    'w_vis_low', 'w_vis_med', 'w_vis_high', 'w_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    for iParam = 1:mp.ncnd
                        k_vis = theta(actparam - 1 + iParam);
                        mp.fulltheta{icnd}.(readparams{iParam}) = k_vis*mp.unifulltheta{unicnd}.(readparams{iParam});
                        mp.fulltheta{icnd}.(readparams{iParam+4}) = k_vis*mp.unifulltheta{unicnd}.(readparams{iParam+4});
                        mp.fulltheta{icnd}.(['k_vis_' num2str(iParam)]) = k_vis; 
                    end
                    k_vest = theta(actparam + mp.ncnd);                    
                    mp.fulltheta{icnd}.sigma_vest = k_vest*mp.unifulltheta{unicnd}.sigma_vest;
                    mp.fulltheta{icnd}.w_vest = k_vest*mp.unifulltheta{unicnd}.w_vest;
                    mp.fulltheta{icnd}.k_vest = k_vest;
                    actparam = actparam + mp.ncnd + 1;
                end
            case 9 % Independently fully proportional to unimodal data, free w (2 params)
                readparams = {'sigma_vis_low', 'sigma_vis_med', 'sigma_vis_high', 'sigma_vest'};
                unicnd = mp.cnd(icnd) - 4; % Look at the same level of noise
                if mp.sharedparams(3) && icnd > 1
                    for iParam = 1:size(readparams, 2)
                        mp.fulltheta{icnd}.(readparams{iParam}) = mp.fulltheta{1}.(readparams{iParam});
                    end
                else
                    for iParam = 1:3
                        mp.fulltheta{icnd}.(readparams{iParam}) = theta(actparam)*mp.unifulltheta{unicnd}.(readparams{iParam});
                    end
                    mp.fulltheta{icnd}.sigma_vest = theta(actparam+1)*mp.unifulltheta{unicnd}.sigma_vest;
                    mp.fulltheta{icnd}.k_vis = theta(actparam);
                    mp.fulltheta{icnd}.k_vest = theta(actparam+1);
                    actparam = actparam + 2;
                end
        end
        
        % Visual internal sensory likelihood
        switch model(4)
            case 1 % Equal to external
                mp.fulltheta{icnd}.sigmalike_vis_low = mp.fulltheta{icnd}.sigma_vis_low;
                mp.fulltheta{icnd}.wlike_vis_low = mp.fulltheta{icnd}.w_vis_low;
                mp.fulltheta{icnd}.sigmalike_vis_med = mp.fulltheta{icnd}.sigma_vis_med;
                mp.fulltheta{icnd}.wlike_vis_med = mp.fulltheta{icnd}.w_vis_med;
                mp.fulltheta{icnd}.sigmalike_vis_high = mp.fulltheta{icnd}.sigma_vis_high;
                mp.fulltheta{icnd}.wlike_vis_high = mp.fulltheta{icnd}.w_vis_high;
        end
        
        % Vestibular internal sensory likelihood
        switch model(5)
            case 1 % Equal to external
                mp.fulltheta{icnd}.sigmalike_vest = mp.fulltheta{icnd}.sigma_vest;
                mp.fulltheta{icnd}.wlike_vest = mp.fulltheta{icnd}.w_vest;
        end
        
        % Visual and vestibular additional noise parameters
        for i = 6:7
            switch model(i)
                case 1 % None                    
                case 2 % Rescaling
                    updateparams(icnd, i, {'id'});
                case 3 % Rescaling and shift
                    updateparams(icnd, i, {'id', 'id'});
                case 5 % Fixed rescaling
            end
        end
        
        % Prior model
        switch model(8)
            case 1 % Gaussian (2-params)
                updateparams(icnd, 8, {'id', 'exp'});
            case 2 % Centered Gaussian (1 params)
                updateparams(icnd, 8, {'exp'});
        end

        % Prior model, part two
        switch model(9)
            case {3,5} % Correlated or discrete free Gaussian (1 params)
                updateparams(icnd, 9, {'exp'});
        end
        
        % Decision making model
        switch model(11)
            case 1 % Standard BDT            
            case 2 % Softmax
                updateparams(icnd, 11, {'exp'});
        end

        % Lapse model
        switch model(13)
            case 1 % No lapse
            case 2 % Prior-matching lapse (1-param)
                updateparams(icnd, 13, {'id'});
        end
        
        % Only for bimodal datasets
        if any(mp.cnd > 4)        
            % Causal inference model (Bimodal only)

            % Model weighting model (Bimodal only)
            switch model(15)
                case {1,6,8,9} % Bayesian posterior
                    updateparams(icnd, 15, {'id'});
                case 2
                    updateparams(icnd, 15, {'id','exp'});
                case 3 % Criterion on x
                    updateparams(icnd, 15, {'exp'});
                case 4 % Soft criterion on x
                    updateparams(icnd, 15, {'exp','exp'});
            end
            
            % Report of unity criterion (Bimodal only)
            switch model(16)
                case 1  % Default (copy localization criterion)
                    mp.fulltheta{icnd}.pcommon_unity = mp.fulltheta{icnd}.pcommon;
                    mp.fulltheta{icnd}.kcommon_unity = mp.fulltheta{icnd}.kcommon;
                    if isfield(mp.fulltheta{icnd},'invgamma_causinf')
                        mp.fulltheta{icnd}.invgamma_causinf_unity = mp.fulltheta{icnd}.invgamma_causinf;
                    end
                case 2  % Separate criterion parameter
                    switch model(15)
                        case {1,2,6,7,8,9} % Bayesian
                            updateparams(icnd, 16, {'id'});
                        case {3,4} % Criterion on x
                            updateparams(icnd, 16, {'exp'});
                    end
                    if isfield(mp.fulltheta{icnd},'invgamma_causinf')
                        mp.fulltheta{icnd}.invgamma_causinf_unity = mp.fulltheta{icnd}.invgamma_causinf;
                    end
                case 3  % Separate criterion and gamma parameter
                    switch model(15)
                        case {1,2,6,7,8,9} % Bayesian
                            updateparams(icnd, 16, {'id','exp'});
                        case {3,4} % Criterion on x
                            updateparams(icnd, 16, {'exp','exp'});
                    end
                case 4 % Random unity
                    updateparams(icnd, 16, {'id','id','id'});
                case 5 % Separate criterion and priorsigmadelta parameter
                    switch model(15)
                        case {1,2,6,7,8,9} % Bayesian
                            updateparams(icnd, 16, {'id','exp'});
                        case {3,4} % Criterion on x
                            updateparams(icnd, 16, {'exp','exp'});
                    end
                case 6 % Forced fusion on localization
                    switch model(15)
                        case {1,2,6,7,8,9} % Bayesian
                            mp.fulltheta{icnd}.pcommon_unity = mp.fulltheta{icnd}.pcommon;
                        case {3,4} % Criterion on x
                            mp.fulltheta{icnd}.kcommon_unity = mp.fulltheta{icnd}.kcommon;
                    end
            end
        end

        
        % mp.fulltheta{1}
    end
    
    % Forced fusion joint fits
    for icnd = 1:mp.ncnd
        if model(16) == 6
            switch model(15)
                case {1,2,6,7,8,9} % Bayesian
                    mp.fulltheta{icnd}.pcommon = 1;
                case {3,4} % Criterion on x
                    mp.fulltheta{icnd}.kcommon = Inf;
            end            
        end
    end
   
    return;
    
    % UPDATEPARAMS Generic parameter update.
    function updateparams(icnd, iParam, oplist)        
        if mp.sharedparams(iParam) && icnd > 1
            for kk = 1:length(mp.paramstree{iParam})
                mp.fulltheta{icnd}.(mp.paramstree{iParam}{kk}) = mp.fulltheta{1}.(mp.paramstree{iParam}{kk});
            end
        else
            for kk = 1:length(mp.paramstree{iParam})
                th = theta(actparam - 1 + kk);
                switch oplist{kk}
                    case 'id' % Identity
                    case 'exp' % Exponential (for log representation)
                        th = exp(th);
                    case 'inv' % Inverse representation
                        th = 1/th - 1;
                end
                mp.fulltheta{icnd}.(mp.paramstree{iParam}{kk}) = th;
            end
            actparam = actparam + length(mp.paramstree{iParam});            
        end
    end
    
end