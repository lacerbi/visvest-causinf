% CUEBMS_BIMODALDATALIKE Calculate (log) likelihood of bimodal dataset X
% under model MODEL and parameter set THETA for binned responses.
%
% X is a unimodal data matrix. Each row is a trial. For a given row, the 
% columns contain data for:
% X(1) Number of trial (unused),
% X(2) Task type (1 Visual localisation, 2 Audio localisation, 3 Categorisation)
% X(3) Audio stimulus position (deg),
% X(4) Visual stimulus position (deg),
% X(5) Response (deg or category).
%
% MODEL is the model class parameter vector.
%
% THETA is the model-dependent parameter vector for the current condition. 
% The values of THETA are:
% THETA(1) sigmazero; THETA(2) szero; 
% THETA(3) alpha (scaling factor); THETA(4) beta (shift factor); 
% THETA(5) kappa (power-law); THETA(6) lambda (lapse rate).
%
% PRIORINFO is the model-depent parameter vector for the prior:
% PRIORINFO(1)  Gaussian mean
% PRIORINFO(2)  Gaussian SD
% PRIORINFO(3)  p_common
% PRIORINFO(4)  k_common
% 
% Example of external usage:
% d = data{1, 3}; model = [6 3 5 2 1]; theta = [0.03 -0.055, 0.06 0.14, 10 0 0.2];
% ParticleCatch_datalike(d.niceData, d.priormix, model, theta, d.priorsinglegauss)
%
function varargout = VestBMS_BimodalLeftRightDatalike_discrete(X,model,theta,priorinfo,bincenters,XGRID,sumover,randomize)

DEBUG = 0;  % Plot some debug graphs

% Program constants
if nargin < 6 || isempty(XGRID) || isnan(XGRID); XGRID = 401; end
if nargin < 7 || isempty(sumover); sumover = 1; end
if nargin < 8 || isempty(randomize); randomize = 0; end

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1e-4;

% Take model parameters
sigmazero_vis = theta(1);
w_vis = theta(2);
sigmalikezero_vis = theta(3);
wlike_vis = theta(4);
alpha_rescaling_vis = 1;

sigmazero_vest = theta(7);
w_vest = theta(8);
sigmalikezero_vest = theta(9);
wlike_vest = theta(10);
alpha_rescaling_vest = 1;

beta_softmax = theta(14);
% This is effectively infinity but solves some numerical instabilities
if (beta_softmax == Inf); beta_softmax = 1e4; end

gamma_causinf = theta(15);
tau_causinf = theta(16);

gamma_causinf_unity = theta(17);
lambda = theta(18);

% Correlated prior
priorsigmadelta = priorinfo(3);

if priorsigmadelta == 0
    error('Discrete prior should have nonzero PRIORSIGMADELTA.');
end

% Model selection parameter
priorc1 = priorinfo(end-3);
kcommon = priorinfo(end-2);
priorc1_unity = priorinfo(end-1);
kcommon_unity = priorinfo(end);

MAXRNG_XMEAS = 180;

% Trials to be computed
do_estimation = ~isempty(X{2});
do_unity = ~isempty(X{3});

% Compute fixed criterion analytically
fixed_criterion_analytic = do_unity && ~do_estimation && model(15) == 3;

% Random unity judgments
if model(16) == 4
    if do_estimation; error('Estimation not supported with random unity judgements.'); end
    if ~do_unity; error('Unity trials must be present with random unity judgement model.'); end

    prmat = [];
    prmat_unity = zeros(numel(bincenters{2}), 2);
    prmat_unity(:,1) = theta(17);
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
    
    [ll,extras] = finalize(prmat,prmat_unity,X,FIXEDLAPSEPDF,nargout > 1,sumover);
    varargout{1} = ll;
    if nargout > 1; varargout{2} = extras; end
    return;
end

% Use distinct criteria for localization vs unity judgement
distinct_criteria = ...
    (priorc1_unity ~= priorc1 || kcommon_unity ~= kcommon) && do_unity;

% Bin centers is a column vector
bincenters_vis = bincenters{1};
bincenters_vest = bincenters{2};

% Compute sensory noise std per trial for vision
if w_vis >= 0; noisemodel_vis = 'A'; else noisemodel_vis = 'C'; end
sigmas_vis = VestBMS_sensoryNoise(noisemodel_vis,bincenters_vis,sigmazero_vis,w_vis);
if isscalar(sigmas_vis); sigmas_vis = sigmas_vis*ones(numel(bincenters_vis),1); end

% Compute sensory noise std per trial for vestibular
if w_vest >= 0; noisemodel_vest = 'A'; else noisemodel_vest = 'C'; end
sigmas_vest = VestBMS_sensoryNoise(noisemodel_vest,bincenters_vest,sigmazero_vest,w_vest);
if isscalar(sigmas_vest); sigmas_vest = sigmas_vest*ones(numel(bincenters_vest),1); end

% Measurements
xrange_vis = zeros(1, XGRID, 1);
xrange_vest = zeros(1, 1, XGRID);
xrange_vis(1, :, 1) = alpha_rescaling_vis*linspace(max(min(bincenters_vis-MAXSD*sigmas_vis),-MAXRNG_XMEAS), min(max(bincenters_vis+MAXSD*sigmas_vis), MAXRNG_XMEAS), XGRID);
xrange_vest (1, 1, :) = alpha_rescaling_vest*linspace(max(min(bincenters_vest-MAXSD*sigmas_vest),-MAXRNG_XMEAS), min(max(bincenters_vest+MAXSD*sigmas_vest), MAXRNG_XMEAS), XGRID);
dx_vis = xrange_vis(1, 2, 1) - xrange_vis(1, 1, 1);
dx_vest = xrange_vest(1, 1, 2) - xrange_vest(1, 1, 1);

% Wrap large noisy measurement around circle
wraparound = MAXRNG_XMEAS >= 180 && ...
        ( min(bincenters_vis-MAXSD*sigmas_vis) <= -180 || max(bincenters_vis+MAXSD*sigmas_vis) >= 180 || ...
        min(bincenters_vest-MAXSD*sigmas_vest) <= -180 || max(bincenters_vest+MAXSD*sigmas_vest) >= 180);

% Add random jitter to XRANGE's to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize
    xrange_vis = xrange_vis + dx_vis*(rand()-0.5);
    xrange_vest = xrange_vest + dx_vest*(rand()-0.5);
end

if nargout > 1 % Save variables for debug or data generation
    extras.xrange_vis = xrange_vis;
    extras.xrange_vest = xrange_vest;
    extras.srange = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type

srange_vis = bincenters_vis(:);
srange_vest = bincenters_vest(:);
idx = srange_vis == srange_vest;
srange_uni = srange_vis(idx);           % C = 1
srange_vis = srange_vis(~idx);          % C = 2, vis
srange_vest = srange_vest(~idx);        % C = 2, vest

% Compute sensory likelihood std for vision
if wlike_vis >= 0; likemodel_vis = 'A'; else likemodel_vis = 'C'; end
sigmasprime_vis = VestBMS_sensoryNoise(likemodel_vis,srange_vis,sigmalikezero_vis,wlike_vis);
sigmasprime_vis_uni = VestBMS_sensoryNoise(likemodel_vis,srange_uni,sigmalikezero_vis,wlike_vis);

% Compute sensory likelihood std for vestibular
if wlike_vest >= 0; likemodel_vest = 'A'; else likemodel_vest = 'C'; end
sigmasprime_vest = VestBMS_sensoryNoise(likemodel_vest,srange_vest,sigmalikezero_vest,wlike_vest);
sigmasprime_vest_uni = VestBMS_sensoryNoise(likemodel_vest,srange_uni,sigmalikezero_vest,wlike_vest);

if fixed_criterion_analytic
    % Fixed-criterion unity judgment
    error('Analytic fixed-criterion model not supported with discrete prior.');
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood for non-Gaussian likelihoods

    if wraparound
        like_vis = bsxfun_normpdf(xrange_vis,srange_vis,sigmasprime_vis) + ...
            bsxfun_normpdf(xrange_vis,srange_vis + 360,sigmasprime_vis) + ...
            bsxfun_normpdf(xrange_vis,srange_vis - 360,sigmasprime_vis);
        like_vest = bsxfun_normpdf(xrange_vest,srange_vest,sigmasprime_vest) + ...
            bsxfun_normpdf(xrange_vest,srange_vest + 360,sigmasprime_vest) + ...
            bsxfun_normpdf(xrange_vest,srange_vest - 360,sigmasprime_vest);
        like_vis_uni = bsxfun_normpdf(xrange_vis,srange_uni,sigmasprime_vis_uni) + ...
            bsxfun_normpdf(xrange_vis,srange_uni + 360,sigmasprime_vis_uni) + ...
            bsxfun_normpdf(xrange_vis,srange_uni - 360,sigmasprime_vis_uni);
        like_vest_uni = bsxfun_normpdf(xrange_vest,srange_uni,sigmasprime_vest_uni) + ...
            bsxfun_normpdf(xrange_vest,srange_uni + 360,sigmasprime_vest_uni) + ...
            bsxfun_normpdf(xrange_vest,srange_uni - 360,sigmasprime_vest_uni);
    else
        like_vis = bsxfun_normpdf(xrange_vis,srange_vis,sigmasprime_vis);
        like_vest = bsxfun_normpdf(xrange_vest,srange_vest,sigmasprime_vest);
        like_vis_uni = bsxfun_normpdf(xrange_vis,srange_uni,sigmasprime_vis_uni);
        like_vest_uni = bsxfun_normpdf(xrange_vest,srange_uni,sigmasprime_vest_uni);
    end
    
    %----------------------------------------------------------------------
    if DEBUG % Plot likelihoods
        subplot(1,2,1); hold off;
        for i = 1:3:numel(srange_vis); plot(xrange_vis(:),like_vis(i,:),'k','LineWidth',1); hold on; end
        subplot(1,2,2); hold off;
        for i = 1:3:numel(srange_vest); plot(xrange_vest(:),like_vest(i,:),'k','LineWidth',1); hold on; end
        stdfig();
        pause
    end
    %----------------------------------------------------------------------

    % Compute C=1 prior, p(s_bar)
    if isfinite(priorinfo(2))
        priorpdf1d = bsxfun_normpdf(srange_uni,priorinfo(1),priorinfo(2));
    else
        priorpdf1d = ones(size(srange_uni));
    end
    priorpdf1d = priorpdf1d/sum(priorpdf1d, 1); % Normalize discrete prior
        
    % Compute unnormalized posterior and rightward posterior (C = 2)
    postpdf_c2_uni = bsxfun(@times, priorpdf1d, like_vest_uni);

    postright_c2 = [];
    % Compute C=2 prior, p(s_vis, s_vest), for eccentric noise
    if isfinite(priorinfo(2))
        priorpdf2d = bsxfun_normpdf(0.5*srange_vest, -0.5*srange_vis + priorinfo(1), priorinfo(2)); 
    else
        priorpdf2d = ones(size(srange_vest));        
    end
    if isfinite(priorsigmadelta)
        priorpdf2d = priorpdf2d .* bsxfun_normpdf(srange_vest, srange_vis, priorsigmadelta);
    end
    priorpdf2d = priorpdf2d/sum(priorpdf2d); % Normalize discrete prior

    %----------------------------------------------------------------------
    if DEBUG % Plot priors
        subplot(1,2,1); hold off;
        plot(srange_uni,priorpdf1d,'k','LineWidth',1);
        text(0.1,0.95,['\sigma = ' num2str(priorinfo(2),'%.1f') ' deg'],'Units','Normalized','FontSize',14);
        subplot(1,2,2); hold off;
        plot(priorpdf2d,'k','LineWidth',1);
        text(0.1,0.95,['\Delta = ' num2str(priorsigmadelta,'%.1f') ' deg'],'Units','Normalized','FontSize',14);
        stdfig();
        pause
    end
    %----------------------------------------------------------------------
    
    likec1 = [];
    likec2 = [];
    if do_estimation
        [postright_c2(1,:,:),likec2(1,:,:)] = VestBMS_c2corrpostandlikec2sum_discrete(priorpdf2d,like_vis,like_vest,srange_vest);
        likec2 = likec2 + realmin; % No volume element, discrete distributions
    end

    % Compute unnormalized posterior and rightward posterior (C = 1)
    if priorc1 > 0
        if do_estimation
            [postright_c1(1,:,:),likec1(1,:,:)] = VestBMS_c1postandlikec1sum_discrete(postpdf_c2_uni, like_vis_uni, srange_uni);
            likec1 = likec1 + realmin;
        else
            postright_c1 = [];
        end
    else
        postright_c1 = 0;
    end

    if nargout > 1 
        extras.postright_c1 = postright_c1;
        extras.postright_c2 = postright_c2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute causal inference weights

    % Compute marginal likelihood, p(x_vis, x_vest|C)

    if model(15) == 1 || model(15) == 2 || model(15) == 6 || model(15) == 7 % (Generalized) Bayesian posterior
        % CASE C=2, Distinct likelihoods, DISCRETE CORRELATED prior
        if isempty(likec2)
            likec2(1,:,:) = VestBMS_likec2corrsum_discrete(priorpdf2d,like_vis,like_vest) + realmin;
        end
        
        % CASE C=1, Unity likelihoods
        if isempty(likec1)
            likec1(1,:,:) = VestBMS_likec1sum_discrete(postpdf_c2_uni,like_vis_uni) + realmin;
        end
    end

    % Compute weight for cue fusion
    switch model(15)
        case {1, 6, 7}  % Model weight is Bayesian posterior p(C=1|x_1, x_2)
            % likec1 = squeeze(likec1);
            w1 = likec1*priorc1./(likec1*priorc1 + likec2*(1-priorc1));
            if distinct_criteria
                w1_unity = likec1*priorc1_unity./(likec1*priorc1_unity + likec2*(1-priorc1_unity));
            end

        case 2  % Generalized Bayesian causal inference
            % likec1 = squeeze(likec1);
            lratio = log(likec1*priorc1) - log(likec2*(1-priorc1));
            w1 = 1./(1 + exp(-gamma_causinf.*lratio));
            if distinct_criteria
                lratio = log(likec1*priorc1_unity) - log(likec2*(1-priorc1_unity));
                w1_unity = 1./(1 + exp(-gamma_causinf_unity.*lratio));
            end

        case 3 % Non-Bayesian, criterion on x
            w1 = abs(bsxfun(@minus, xrange_vis, xrange_vest)) < kcommon;
            if wraparound
                w1 = w1 | (abs(bsxfun(@minus, xrange_vis + 360, xrange_vest)) < kcommon) ...
                    | (abs(bsxfun(@minus, xrange_vis - 360, xrange_vest)) < kcommon);
            end
            if distinct_criteria
                w1_unity = abs(bsxfun(@minus, xrange_vis, xrange_vest)) < kcommon_unity;
                if wraparound
                    w1_unity = w1_unity | (abs(bsxfun(@minus, xrange_vis + 360, xrange_vest)) < kcommon_unity) ...
                        | (abs(bsxfun(@minus, xrange_vis - 360, xrange_vest)) < kcommon_unity);                    
                end
            end

        case 4 % Soft criterion on x
            w1 = 1./(1 + exp(tau_causinf*(abs(bsxfun(@minus, xrange_vis, xrange_vest)) - kcommon)));
            if distinct_criteria
                w1_unity = 1./(1 + exp(tau_causinf*(abs(bsxfun(@minus, xrange_vis, xrange_vest)) - kcommon_unity)));                
            end

        case 5 % Forced fusion
            w1 = 1;
            
    end

    % NaNs can emerge as 0/0 - assume that the response becomes random
    w1(isnan(w1)) = 0.5;

    % Bisensory estimation
    if do_estimation
        if model(15) == 6   % Posterior model probability matching
            % Compute the probability of each causal structure separately
            postright_c1 = 1./(1 + ((1-postright_c1)./postright_c1).^beta_softmax);
            postright_c1(isnan(postright_c1)) = 0.5;
            postright_c2 = 1./(1 + ((1-postright_c2)./postright_c2).^beta_softmax);
            postright_c2(isnan(postright_c2)) = 0.5;
            % Combine with probability matching
            prright = bsxfun(@times, w1, postright_c1) + bsxfun(@times, 1 - w1, postright_c2);
        else
            % Compute posterior probability of rightward motion    
            postright = bsxfun(@plus, bsxfun(@times, w1, postright_c1), bsxfun(@times, 1-w1, postright_c2));
            
            % Probability of rightward response            
            prright = 1./(1 + ((1-postright)./postright).^beta_softmax);
            prright(isnan(prright)) = 0.5;
        end
    end

    % Bisensory unity judgement
    if do_unity && ~distinct_criteria
        if beta_softmax == 1
            w1_unity = w1;
        elseif beta_softmax == Inf
            w1_unity = zeros(size(w1));
            w1_unity(w1 > 0.5) = 1;
            w1_unity(w1 == 0.5) = 0.5;
        else
            w1_unity = 1./(1 + ((1-w1)./w1).^beta_softmax);
        end
    elseif do_unity
        w1_unity = double(w1_unity);    % Convert from logical to double
        if beta_softmax ~= 1
            w1_unity = 1./(1 + ((1-w1_unity)./w1_unity).^beta_softmax);            
        end
    end

    if distinct_criteria; w1_unity(isnan(w1_unity)) = 0.5; end

    if nargout > 1 % Save variables for debug or data generation
        extras.w1 = w1;
        if distinct_criteria; extras.w1_unity = w1_unity; end
        % extras.postpdfc1 = postpdfc1;
    end
    
    % Clean up memory
    clear postright w1 postpdf_c2_uni postright_c1 postright_c2 ...
        likec1 likec2 likec2_vis postpdfc1 lratio;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marginalize over noisy measurements

if ~fixed_criterion_analytic
    xpdf_vis = bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis,alpha_rescaling_vis*sigmas_vis);
    if wraparound
        xpdf_vis = xpdf_vis + ...
            bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis + 360,alpha_rescaling_vis*sigmas_vis) ...
            + bsxfun_normpdf(xrange_vis, alpha_rescaling_vis*bincenters_vis - 360,alpha_rescaling_vis*sigmas_vis);        
    end
    xpdf_vis = bsxfun(@rdivide, xpdf_vis, qtrapz(xpdf_vis, 2)); % Not multiplying by volume element
    do_symmetrize = 0;
    if do_symmetrize
        mid = ceil(0.5*size(xpdf_vis,1));
        xpdf_vis(mid,:) = 0.5*(xpdf_vis(mid,:) + fliplr(xpdf_vis(mid,:)));
    end
    
    xpdf_vest = bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest,alpha_rescaling_vest*sigmas_vest);    
    if wraparound
        xpdf_vest = xpdf_vest + ...
            bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest + 360,alpha_rescaling_vest*sigmas_vest) ...
            + bsxfun_normpdf(xrange_vest, alpha_rescaling_vest*bincenters_vest - 360,alpha_rescaling_vest*sigmas_vest);        
    end    
    xpdf_vest = bsxfun(@rdivide, xpdf_vest, qtrapz(xpdf_vest, 3));  % Not multiplying by volume element
    if do_symmetrize
        xpdf_vest(mid,1,:) = 0.5*(xpdf_vest(mid,1,:) + xpdf_vest(mid,1,end:-1:1));
    end
    
    %----------------------------------------------------------------------
    if DEBUG % Plot noisy measurements distributions
        subplot(1,2,1); hold off;
        for i = 1:3:numel(bincenters_vis); plot(xrange_vis(:),xpdf_vis(i,:),'k','LineWidth',1); hold on; end
        subplot(1,2,2); hold off;
        for i = 1:3:numel(bincenters_vest); plot(xrange_vest(:),xpdf_vest(i,:),'k','LineWidth',1); hold on; end
        stdfig();
        pause
    end
    %----------------------------------------------------------------------    

    if do_estimation
        prmat = zeros(numel(bincenters_vest), 2);
        prmat(:,2) = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,prright);   % Not multiplying by volume element (xpdfs did not)
        prmat(:,1) = 1 - prmat(:,2);
    else
        prmat = [];
    end

    if do_unity
        prmat_unity = zeros(numel(bincenters_vest), 2);
        prmat_unity(:,1) = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,w1_unity);    % Not multiplying by volume element (xpdfs did not)
        prmat_unity(:,2) = 1 - prmat_unity(:,1);
    else
        prmat_unity = [];
    end
end

% Fix probabilities
prmat = min(max(prmat,0),1);
prmat_unity = min(max(prmat_unity,0),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize log likelihood

prmat = lambda/2 + (1-lambda)*prmat;
prmat_unity = lambda/2 + (1-lambda)*prmat_unity;

[ll,extras] = finalize(prmat,prmat_unity,X,FIXEDLAPSEPDF,nargout > 1,sumover);
varargout{1} = ll;
if nargout > 1; varargout{2} = extras; end

end

%--------------------------------------------------------------------------
function [ll,extras] = finalize(prmat,prmat_unity,X,epsilon,extrasflag,sumoverflag)
%FINALIZE Finalize log likelihood

prmat = 0.5*epsilon + (1-epsilon)*prmat;
prmat_unity = 0.5*epsilon + (1-epsilon)*prmat_unity;

if extrasflag
    extras.responsepdf = prmat;
    extras.responsepdf_unity = prmat_unity;
else
    extras = [];
end

prmat = [prmat(:); prmat_unity(:)];
xx = [X{1}(:); X{2}(:); X{3}(:)];

if sumoverflag
    ll = sum(xx.*log(prmat));
else
    ll = loglikmat2vec(log(prmat),xx);
end

end