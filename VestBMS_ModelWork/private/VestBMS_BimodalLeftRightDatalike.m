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
function varargout = VestBMS_BimodalLeftRightDatalike(X,model,theta,priorinfo,bincenters,maxranges,XGRID,SSCALE,sumover,randomize)

DEBUG = 0;  % Plot some debug graphs

% Program constants
if nargin < 7 || isempty(XGRID) || isnan(XGRID); XGRID = 401; end
if nargin < 8 || isempty(SSCALE); SSCALE = 8; end
if nargin < 9 || isempty(sumover); sumover = 1; end
if nargin < 10 || isempty(randomize); randomize = 0; end

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
% beta_rescaling_vis = 0;

sigmazero_vest = theta(7);
w_vest = theta(8);
sigmalikezero_vest = theta(9);
wlike_vest = theta(10);
alpha_rescaling_vest = 1;
% beta_rescaling_vest = 0;

beta_softmax = theta(14);
% This is effectively infinity but solves some numerical instabilities
if (beta_softmax == Inf); beta_softmax = 1e4; end

gamma_causinf = theta(15);
tau_causinf = theta(16);

gamma_causinf_unity = theta(17);
lambda = theta(18);

% Correlated prior
if numel(priorinfo) < 7; priorsigmadelta = 0; else priorsigmadelta = priorinfo(3); end

% Model selection parameter
priorc1 = priorinfo(end-3);
kcommon = priorinfo(end-2);
priorc1_unity = priorinfo(end-1);
kcommon_unity = priorinfo(end);

MAXRNG = maxranges(1);
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
    
    [ll,extras] = finalize(prmat,prmat_unity,X,FIXEDLAPSEPDF,nargout > 1,sumover(1));
    varargout{1} = ll;
    if nargout > 1; varargout{2} = extras; end
    return;
end

% Use distinct criteria for localization vs unity judgement
distinct_criteria = ...
    (priorc1_unity ~= priorc1 || kcommon_unity ~= kcommon) && do_unity;

% Is the internal noise Gaussian?
gaussianflag = (wlike_vis == 0 && wlike_vest == 0);

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

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
ds = diff(srange(1:2));

if nargout > 1 % Save variables for debug or data generation
    extras.xrange_vis = xrange_vis;
    extras.xrange_vest = xrange_vest;
    extras.srange = srange;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type

if gaussianflag
    likerange_vis = xrange_vis; 
    likerange_vest = xrange_vest; 
else
    likerange_vis = srange; 
    likerange_vest = srange; 
end

% Compute sensory likelihood std for vision
if wlike_vis >= 0; likemodel_vis = 'A'; else likemodel_vis = 'C'; end
sigmasprime_vis = VestBMS_sensoryNoise(likemodel_vis,likerange_vis,sigmalikezero_vis,wlike_vis);

% Compute sensory likelihood std for vestibular
if wlike_vest >= 0; likemodel_vest = 'A'; else likemodel_vest = 'C'; end
sigmasprime_vest = VestBMS_sensoryNoise(likemodel_vest,likerange_vest,sigmalikezero_vest,wlike_vest);

if fixed_criterion_analytic
    % Fixed-criterion unity judgment
    
    mu_diff = bincenters_vest - bincenters_vis;
    sigma_diff = sqrt(sigmas_vis.^2 + sigmas_vest.^2);
    
    prmat = [];
    prmat_unity = zeros(numel(bincenters_vest), 2);    
    prmat_unity(:,1) = bsxfun_normcdf(kcommon,mu_diff,sigma_diff) - bsxfun_normcdf(-kcommon,mu_diff,sigma_diff);
    if wraparound
        prmat_unity(:,1) = prmat_unity(:,1) + bsxfun_normcdf(kcommon,mu_diff-360,sigma_diff) - bsxfun_normcdf(-kcommon,mu_diff-360,sigma_diff);
        prmat_unity(:,1) = prmat_unity(:,1) + bsxfun_normcdf(kcommon,mu_diff+360,sigma_diff) - bsxfun_normcdf(-kcommon,mu_diff+360,sigma_diff);
    end    
    prmat_unity(:,2) = 1 - prmat_unity(:,1);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood for non-Gaussian likelihoods

    if do_estimation || ~gaussianflag
        if wraparound
            like_vis = bsxfun_normpdf(xrange_vis,srange,sigmasprime_vis) + ...
                bsxfun_normpdf(xrange_vis,srange + 360,sigmasprime_vis) + ...
                bsxfun_normpdf(xrange_vis,srange - 360,sigmasprime_vis);
            like_vest = bsxfun_normpdf(xrange_vest,srange,sigmasprime_vest) + ...
                bsxfun_normpdf(xrange_vest,srange + 360,sigmasprime_vest) + ...
                bsxfun_normpdf(xrange_vest,srange - 360,sigmasprime_vest);
        else
            like_vis = bsxfun_normpdf(xrange_vis,srange,sigmasprime_vis);
            like_vest = bsxfun_normpdf(xrange_vest,srange,sigmasprime_vest);
        end
    
        % Compute UNCORRELATED prior, p(s)
        priorpdf1d = bsxfun_normpdf(srange,priorinfo(1),priorinfo(2));
        priorpdf1d = priorpdf1d/(qtrapz(priorpdf1d, 1)*ds); % Normalize prior

        % Compute unnormalized posterior and rightward posterior (C = 2)
        postpdf_c2 = bsxfun(@times, priorpdf1d, like_vest);
        
        postright_c2 = [];
        if priorsigmadelta > 0 && ~gaussianflag
            
            adaptiveintegration = 0;
            
            % Compute CORRELATED prior, p(s_vis, s_vest), for eccentric noise
            if adaptiveintegration
                srange_vis_adapt = linspace(-MAXRNG, MAXRNG, 1e3)';
                srange_vest_adapt = srange_vis_adapt(:)';
                
                % Compute sensory likelihood std for vision
                sigmasprime_vis_adapt = VestBMS_sensoryNoise(likemodel_vis,srange_vis_adapt,sigmalikezero_vis,wlike_vis);

                % Compute sensory likelihood std for vestibular
                sigmasprime_vest_adapt = VestBMS_sensoryNoise(likemodel_vest,srange_vest_adapt,sigmalikezero_vest,wlike_vest);
                                
                if wraparound
                    like_vis_adapt = bsxfun_normpdf(xrange_vis,srange_vis_adapt,sigmasprime_vis_adapt) + ...
                        bsxfun_normpdf(xrange_vis,srange_vis_adapt + 360,sigmasprime_vis_adapt) + ...
                        bsxfun_normpdf(xrange_vis,srange_vis_adapt - 360,sigmasprime_vis_adapt);
                    like_vest_adapt(:,1,:) = bsxfun_normpdf(xrange_vest,srange_vest_adapt,sigmasprime_vest_adapt) + ...
                        bsxfun_normpdf(xrange_vest,srange_vest_adapt + 360,sigmasprime_vest_adapt) + ...
                        bsxfun_normpdf(xrange_vest,srange_vest_adapt - 360,sigmasprime_vest_adapt);
                else
                    like_vis_adapt = bsxfun_normpdf(xrange_vis,srange_vis_adapt,sigmasprime_vis_adapt);
                    like_vest_adapt(:,1,:) = bsxfun_normpdf(xrange_vest,srange_vest_adapt,sigmasprime_vest_adapt);
                end
                priorpdf2d = bsxfun_normpdf(0.5*srange_vest_adapt, -0.5*srange_vis_adapt + priorinfo(1), priorinfo(2)) .* ...
                    bsxfun_normpdf(srange_vest_adapt, srange_vis_adapt, priorsigmadelta);
                priorpdf2d = priorpdf2d/(qtrapz(qtrapz(priorpdf2d,1))*ds*ds); % Normalize prior
            else
                srange_vis = srange(:);
                srange_vest = srange(:)';
                priorpdf2d = bsxfun_normpdf(0.5*srange_vest, -0.5*srange_vis + priorinfo(1), priorinfo(2)) .* ...
                    bsxfun_normpdf(srange_vest, srange_vis, priorsigmadelta);
                priorpdf2d = priorpdf2d/(qtrapz(qtrapz(priorpdf2d,1))*ds*ds); % Normalize prior
            end
                        
            if do_estimation
                if adaptiveintegration
                    tic
                    [postright_c2(1,:,:),likec2(1,:,:),err,fevals] = VestBMS_c2corrpostandlikec2adapt_mex2(priorpdf2d,like_vis_adapt,like_vest_adapt,2e-3);
                    toc
                    if 1
                        tic
                        [postright_c2_2(1,:,:),likec2_2(1,:,:),err,fevals] = VestBMS_c2corrpostandlikec2adapt_mex2(priorpdf2d,like_vis_adapt,like_vest_adapt,1e-5);
                        toc
                        srange_vis = srange(:);
                        srange_vest = srange(:)';
                        priorpdf2d = bsxfun_normpdf(0.5*srange_vest, -0.5*srange_vis + priorinfo(1), priorinfo(2)) .* ...
                            bsxfun_normpdf(srange_vest, srange_vis, priorsigmadelta);
                        priorpdf2d = priorpdf2d/(qtrapz(qtrapz(priorpdf2d,1))*ds*ds); % Normalize prior
                        tic
                        [postright_c2_3(1,:,:),likec2_3(1,:,:)] = VestBMS_c2corrpostandlikec2qtrapz(priorpdf2d,like_vis,like_vest);
                        toc
                        pause
                    end
                else
                    [postright_c2(1,:,:),likec2(1,:,:)] = VestBMS_c2corrpostandlikec2qtrapz(priorpdf2d,like_vis,like_vest);
                end
                likec2 = likec2*ds*ds + realmin;
            else
                likec2 = [];
            end
            
        elseif do_estimation
            if priorc1 < 1
                if priorsigmadelta == 0
                    postright_c2(1,:,:) = VestBMS_PostRight(postpdf_c2);
                else % CORRELATED prior with constant noise
                    T = priorsigmadelta^2 * bsxfun(@minus, xrange_vest .* sigmasprime_vis(1)^2, xrange_vis .* sigmasprime_vest(1)^2) ...
                        + 4*priorinfo(2)^2 * bsxfun(@plus, xrange_vest * (priorsigmadelta^2 + sigmasprime_vis(1)^2), xrange_vis * sigmasprime_vest(1)^2);
                    S2 = (priorsigmadelta^2 * sigmasprime_vis(1)^2 + 4 * priorinfo(2)^2 * (priorsigmadelta^2 + sigmasprime_vis(1)^2)) ...
                        * sigmasprime_vest(1)^2 * (4 * sigmasprime_vis(1)^2 * sigmasprime_vest(1)^2 + ...
                        priorsigmadelta^2 * (sigmasprime_vis(1)^2 * sigmasprime_vest(1)^2) + 4 * priorinfo(2)^2 * (priorsigmadelta^2 + sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2));
                    postright_c2 = 0.5 * ( 1 + erf(T / sqrt(2*S2)) );
                end
            else
                postright_c2 = 0;
            end
        end
        
        likec1 = [];

        % Compute unnormalized posterior and rightward posterior (C = 1)
        if priorc1 > 0
            if do_estimation
                if gaussianflag
                    mutilde = bsxfun(@plus, xrange_vest.*sigmasprime_vis(1)^2, xrange_vis.*sigmasprime_vest(1)^2)./(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2);
                    sigma2tilde = sigmasprime_vis(1)^2.*sigmasprime_vest(1)^2./(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2);
                    mucdf = mutilde.*priorinfo(2)^2./(sigma2tilde + priorinfo(2)^2);
                    sigmacdf = sqrt(sigma2tilde./(sigma2tilde + priorinfo(2)^2))*priorinfo(2);
                    
                    intright = (bsxfun_normcdf(MAXRNG, mucdf, sigmacdf) - bsxfun_normcdf(0, mucdf, sigmacdf));
                    if wraparound
                        intright = intright + ...
                            (bsxfun_normcdf(MAXRNG + 360, mucdf, sigmacdf) - bsxfun_normcdf(360, mucdf, sigmacdf)) + ...
                            (bsxfun_normcdf(MAXRNG - 360, mucdf, sigmacdf) - bsxfun_normcdf(-360, mucdf, sigmacdf));
                    end
                    %likeright = intright .* bsxfun_normpdf(xrange_vest,xrange_vis,sqrt(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2)) .* ...
                    %    bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde + priorinfo(2)^2)) + realmin;                
                    
                    intleft = (bsxfun_normcdf(0, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG, mucdf, sigmacdf));
                    if wraparound
                        intleft = intleft + ...
                            (bsxfun_normcdf(360, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG + 360, mucdf, sigmacdf)) + ...
                            (bsxfun_normcdf(-360, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG - 360, mucdf, sigmacdf));
                    end            
                    %likeleft = intleft .* bsxfun_normpdf(xrange_vest,xrange_vis,sqrt(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2)) .* ...
                    %    bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde + priorinfo(2)^2)) + realmin;
                    
                    postright_c1(1,:,:) = intright./(intleft + intright);
                    
                else
                    [postright_c1(1,:,:),likec1(1,:,:)] = VestBMS_c1postandlikec1qtrapz(postpdf_c2, like_vis);
                    likec1 = likec1*ds + realmin;   % ADDED DS!
                    % likec1 = likec1 + realmin;
                end
            else
                postright_c1 = [];
            end
        else
            postright_c1 = 0;
        end
    else
        postright_c1 = [];
        postright_c2 = [];
    end

    if nargout > 1 
        extras.postright_c1 = postright_c1; 
        extras.postright_c2 = postright_c2;    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute causal inference weights

    % Compute marginal likelihood, p(x_vis, x_vest|C)

    if model(15) == 1 || model(15) == 2 || model(15) == 6 || model(15) == 7 || model(15) == 8 % (Generalized) Bayesian posterior
        % CASE C=2, Independent likelihoods
        if gaussianflag
            if priorsigmadelta == 0
                muc2_vis = xrange_vis.*priorinfo(2)^2/(sigmasprime_vis(1)^2 + priorinfo(2)^2);
                sigmac2_vis = sigmasprime_vis(1)*priorinfo(2)/sqrt(sigmasprime_vis(1)^2 + priorinfo(2)^2);
                muc2_vest = xrange_vest.*priorinfo(2)^2/(sigmasprime_vest(1)^2 + priorinfo(2)^2);
                sigmac2_vest = sigmasprime_vest(1)*priorinfo(2)/sqrt(sigmasprime_vest(1)^2 + priorinfo(2)^2);            
                int_vis = (bsxfun_normcdf(MAXRNG,muc2_vis,sigmac2_vis) - bsxfun_normcdf(-MAXRNG,muc2_vis,sigmac2_vis));
                int_vest = (bsxfun_normcdf(MAXRNG,muc2_vest,sigmac2_vest) - bsxfun_normcdf(-MAXRNG,muc2_vest,sigmac2_vest));
                if wraparound
                    int_vis = int_vis + (bsxfun_normcdf(MAXRNG + 360,muc2_vis,sigmac2_vis) - bsxfun_normcdf(-MAXRNG + 360,muc2_vis,sigmac2_vis)) ...
                        + (bsxfun_normcdf(MAXRNG - 360,muc2_vis,sigmac2_vis) - bsxfun_normcdf(-MAXRNG - 360,muc2_vis,sigmac2_vis));
                    int_vest = int_vest + (bsxfun_normcdf(MAXRNG + 360,muc2_vest,sigmac2_vest) - bsxfun_normcdf(-MAXRNG + 360,muc2_vest,sigmac2_vest)) ...
                        + (bsxfun_normcdf(MAXRNG - 360,muc2_vest,sigmac2_vest) - bsxfun_normcdf(-MAXRNG - 360,muc2_vest,sigmac2_vest));                
                end
                int_vis = int_vis .* bsxfun_normpdf(xrange_vis,0,sqrt(sigmasprime_vis(1)^2 + priorinfo(2)^2));
                int_vest = int_vest .* bsxfun_normpdf(xrange_vest,0,sqrt(sigmasprime_vest(1)^2 + priorinfo(2)^2));
                likec2 = bsxfun(@times, int_vis, int_vest) + realmin;
            else
                sigma2star = 4*sigmasprime_vis(1)^2*sigmasprime_vest(1)^2 + ...
                    priorsigmadelta^2*(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2) + ...
                    4* priorinfo(2)^2*(priorsigmadelta^2 + sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2);
                z2 = 4* priorinfo(2)^2 * bsxfun(@minus, xrange_vis, xrange_vest).^2 + ...
                    bsxfun(@plus, ...
                        bsxfun(@plus, 2* bsxfun(@times, xrange_vis, xrange_vest) * priorsigmadelta^2, ...
                        xrange_vest.^2 * (priorsigmadelta^2 + 4 * sigmasprime_vis(1)^2)), ...
                    xrange_vis.^2 * (priorsigmadelta^2 + 4 * sigmasprime_vest(1)^2));
                likec2 = exp(-0.5*z2./sigma2star)./(pi*sqrt(sigma2star));
            end
        elseif priorsigmadelta > 0
            % CORRELATED prior
            if isempty(likec2)
                likec2(1,:,:) = VestBMS_likec2corrqtrapz(priorpdf2d,like_vis,like_vest)*ds*ds + realmin;
            end
        else
            % UNCORRELATED prior
            likec2_vis = bsxfun(@times, priorpdf1d, like_vis);
            likec2 = (bsxfun(@times, qtrapz(likec2_vis, 1)*ds, qtrapz(postpdf_c2, 1)*ds)) + realmin;
        end

        % postpdfc1 = bsxfun(@rdivide, postpdfc1, likec1);
        
        % CASE C=1, Likelihoods are not independent
        if gaussianflag
            if priorinfo(1) ~= 0; error('Prior bias not supported yet.'); end
            mutilde = bsxfun(@plus, xrange_vest.*sigmasprime_vis(1)^2, xrange_vis.*sigmasprime_vest(1)^2)./(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2);
            sigma2tilde = sigmasprime_vis(1)^2.*sigmasprime_vest(1)^2./(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2);
            mucdf = mutilde.*priorinfo(2)^2./(sigma2tilde + priorinfo(2)^2);
            sigmacdf = sqrt(sigma2tilde./(sigma2tilde + priorinfo(2)^2))*priorinfo(2);
            intc1 = (bsxfun_normcdf(MAXRNG, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG, mucdf, sigmacdf));
            if wraparound
                intc1 = intc1 + ...
                    (bsxfun_normcdf(MAXRNG + 360, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG + 360, mucdf, sigmacdf)) + ...
                    (bsxfun_normcdf(MAXRNG - 360, mucdf, sigmacdf) - bsxfun_normcdf(-MAXRNG - 360, mucdf, sigmacdf));
            end            
            likec1 = intc1 .* bsxfun_normpdf(xrange_vest,xrange_vis,sqrt(sigmasprime_vis(1)^2 + sigmasprime_vest(1)^2)) .* ...
                bsxfun_normpdf(mutilde,0,sqrt(sigma2tilde + priorinfo(2)^2)) + realmin;                
        else
            if isempty(likec1)
                likec1(1,:,:) = VestBMS_likec1qtrapz(postpdf_c2,like_vis)*ds + realmin; % ADDED DS!
                %likec1(1,:,:) = VestBMS_likec1qtrapz(postpdf_c2,like_vis) + realmin;
            end
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
            
        case 8 % Posterior model selection
            w1 = double(likec1*priorc1 >= likec2*(1-priorc1));
            if distinct_criteria
                w1_unity = likec1*priorc1_unity./(likec1*priorc1_unity + likec2*(1-priorc1_unity));
            end
            
        case 9 % Probabilistic fusion
            w1 = priorc1;
    end

    % NaNs can emerge as 0/0 - assume that the response becomes random
    w1(isnan(w1)) = 0.5;

    % Bisensory estimation
    if do_estimation
        if model(15) == 6 || model(15) == 9  % Posterior model probability matching or probabilistic fusion
            % Compute the probability of each causal structure separately
            % postright_c1 = 1./(1 + ((1-postright_c1)./postright_c1).^beta_softmax);
            % postright_c1 = 1./(1 + ((1-postright_c1)./postright_c1));
            postright_c1(isnan(postright_c1)) = 0.5;
            % postright_c2 = 1./(1 + ((1-postright_c2)./postright_c2).^beta_softmax);
            % postright_c2 = 1./(1 + ((1-postright_c2)./postright_c2));
            postright_c2(isnan(postright_c2)) = 0.5;
            % Combine with probability matching
            prright = bsxfun(@plus, bsxfun(@times, w1, postright_c1), bsxfun(@times, 1 - w1, postright_c2));
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
        % w1_unity = w1;
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
    
    %----------------------------------------------------------------------
    % Plot decision rule for unity judgments (function of x_vest, x_vis)
    if DEBUG && do_unity
        clf;
        surf(xrange_vis(:), xrange_vest(:), squeeze(w1_unity),'LineStyle','none'); hold on;
        xlabel('$x_\mathrm{vis}$ (deg)','Interpreter','LaTeX','FontSize',14); 
        ylabel('$x_\mathrm{vest}$ (deg)','Interpreter','LaTeX','FontSize',14);
        axis square;
        xylims = [min(xrange_vest(1),xrange_vis(1)),max(xrange_vest(end),xrange_vis(end))];
        xlim(xylims);
        ylim(xylims);
        plot3(xylims,xylims,[200 200],'--k','LineWidth',1);
        view([0 90]);
        hold off;
        title(['$\sigma^0_\mathrm{vis} = ' num2str(sigmazero_vis,'%.2f') '$ deg, $w_\mathrm{vis} = ' num2str(abs(w_vis),'%.2f') '$' ...
            ', $\sigma^0_\mathrm{vest} = ' num2str(sigmazero_vest,'%.2f') '$ deg, $w_\mathrm{vest} = ' num2str(abs(w_vest),'%.2f') '$'], ...
            'FontSize',14,'Interpreter','LaTeX');
        pause;
    end
    %----------------------------------------------------------------------    
    
    % Clean up memory
    clear postright w1 postpdf_c2 postright_c1 postright_c2 ...
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

[ll,extras] = finalize(prmat,prmat_unity,X,FIXEDLAPSEPDF,nargout > 1,sumover(1));
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