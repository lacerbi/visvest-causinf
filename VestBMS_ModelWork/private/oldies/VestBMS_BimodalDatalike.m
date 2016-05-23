% CUEBMS_BIMODALDATALIKE Calculate log likelihood of bimodal dataset X
% under model MODEL and parameter set THETA.
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
% PRIORINFO (1) mean of the first Gaussian
% PRIORINFO (2) SD of the first Gaussian
% PRIORINFO (3) mixing weight of the first Gaussian
% PRIORINFO (4) mean of the second Gaussian
% PRIORINFO (5) SD of the second Gaussian
% 
% Example of external usage:
% d = data{1, 3}; model = [6 3 5 2 1]; theta = [0.03 -0.055, 0.06 0.14, 10 0 0.2];
% ParticleCatch_datalike(d.niceData, d.priormix, model, theta, d.priorsinglegauss)
%
function varargout = CueBMS_BimodalDatalike(X,model,theta,priorinfo,XGRID,SSCALE,sumover,MAXRNG,randomize)

% Program constants
if nargin < 5 || isempty(XGRID) || isnan(XGRID); XGRID = 201; end
if nargin < 6 || isempty(SSCALE); SSCALE = 8; end
if nargin < 7 || isempty(sumover); sumover = 1; end
if nargin < 8 || isempty(MAXRNG); MAXRNG = 45; end % Screen bounds (in deg)
if nargin < 9 || isempty(randomize); randomize = 0; end

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;
if XGRID < 101; MAXSD = 3; end

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1.5e-6;

% This is almost zero
NUMZERO = 1e-80;

% Take model parameters
sigmazero_vis = theta(1);
szero_vis = theta(2);
sigmalikezero_vis = theta(3);
slikezero_vis = theta(4);
alpha_vis = theta(5);
beta_vis = theta(6);

sigmazero_aud = theta(7);
szero_aud = theta(8);
sigmalikezero_aud = theta(9);
slikezero_aud = theta(10);
alpha_aud = theta(11);
beta_aud = theta(12);

sigmaloss = theta(13);
kappa = theta(14);
lambda = theta(15);

sigmar = theta(19); % Motor/placement noise
wr = theta(22); % Weber's fraction for motor/placement noise
priorc1cat = theta(20); % P(C=1) for categorization trials
lambdacat = theta(21); % Categorisation trials lapse
if isnan(lambdacat); lambdacat = lambda; end

rho = 0.5; % Criterion for model selection
softmaxbeta = 400; % Softmax parameter for model selection

% Model selection parameter (meaning depends on model)
priorc1 = priorinfo(end-1);
kcommon = priorinfo(end);

% Is the internal noise Gaussian?
if ( ( (slikezero_vis == 0 || isinf(slikezero_vis)) && ...
        (slikezero_aud == 0 || isinf(slikezero_aud)) ) ...
        || ((model(4) == 6 || model(4) == 8)  && model(5) == 6)) ...
                            && priorinfo(6) == 0
    gaussianflag = 1;
else
    gaussianflag = 0;
    if model(4) == 6 || model(5) == 6; error('Unsupported measurement-based likelihood.'); end
end

% Number of trials in this condition
nTrials = size(X, 1);

% Compute sensory noise std per trial for vision
if isinf(szero_vis) || szero_vis == 0
    sigmas_vis = sigmazero_vis*ones(nTrials, 1);
elseif szero_vis > 0
    sigmas_vis = sigmazero_vis*sqrt(1 + (X(:, 4)/szero_vis).^2);    
else % Negative szero for cosine noise formula
    sigmas_vis = sigmazero_vis.*sqrt(1 + 2*(90/pi)^2*(1 - cos(X(:,4)*pi/90))/szero_vis^2);
    %sigmas_vis = sigmazero_vis*(1 + abs(X(:, 4)/szero_vis));
end

% Compute sensory noise std per trial for audition
if isinf(szero_aud) || szero_aud == 0
    sigmas_aud = sigmazero_aud*ones(nTrials, 1);
elseif szero_aud > 0
    sigmas_aud = sigmazero_aud*sqrt(1 + (X(:, 3)/szero_aud).^2);    
else % Negative szero for cosine noise formula (horrible way, I know)
    sigmas_aud = sigmazero_aud.*sqrt(1 + 2*(90/pi)^2*(1 - cos(X(:,3)*pi/90))/szero_aud^2);
%    sigmas_aud = sigmazero_aud*(1 + abs(X(:, 3)/szero_aud));        
end

% Measurements
xrange_vis = zeros(1, XGRID, 1);
xrange_aud = zeros(1, 1, XGRID);
xrange_vis(1, :, 1) = alpha_vis*linspace(max(min(X(:,4)-MAXSD*sigmas_vis),-MAXRNG), min(max(X(:,4)+MAXSD*sigmas_vis), MAXRNG), XGRID);
xrange_aud (1, 1, :) = alpha_aud*linspace(max(min(X(:,3)-MAXSD*sigmas_aud),-MAXRNG), min(max(X(:,3)+MAXSD*sigmas_aud), MAXRNG), XGRID);
dx_vis = xrange_vis(1, 2, 1) - xrange_vis(1, 1, 1);
dx_aud = xrange_aud(1, 1, 2) - xrange_aud(1, 1, 1);

softmaxbeta_x = 0.5/min([dx_vis,dx_aud]);     % Softmax beta when applied to x

% Add random jitter to XRANGE's to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize
    xrange_vis = xrange_vis + dx_vis*(rand()-0.5);
    xrange_aud = xrange_aud + dx_aud*(rand()-0.5);
end

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
ds = diff(srange(1:2));
responsepdf = NaN(nTrials, 1);

if nargout > 1 % Save variables for debug or data generation
    extras.xrange_vis = xrange_vis;
    extras.xrange_aud = xrange_aud;
    extras.srange = srange;
end

if ~gaussianflag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood for non-Gaussian likelihoods
    
    % Compute sensory likelihood std for vision
    if isinf(slikezero_vis) || slikezero_vis == 0
        sigmasprime_vis = sigmalikezero_vis;
    elseif slikezero_vis > 0
        sigmasprime_vis = sigmalikezero_vis*sqrt(1 + (srange/slikezero_vis).^2);    
    else % Negative szero for cosine noise formula
        sigmasprime_vis = sigmalikezero_vis.*sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90))/slikezero_vis^2);
        % sigmasprime_vis = sigmalikezero_vis*(1 + abs(srange/slikezero_vis));        
    end

    % Compute sensory likelihood std for audition
    if isinf(slikezero_aud) || slikezero_aud == 0
        sigmasprime_aud = sigmalikezero_aud;
    elseif slikezero_aud > 0
        sigmasprime_aud = sigmalikezero_aud*sqrt(1 + (srange/slikezero_aud).^2);    
    else % Negative szero for cosine noise formula
        sigmasprime_aud = sigmalikezero_aud.*sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90))/slikezero_aud^2);
        % sigmasprime_aud = sigmalikezero_aud*(1 + abs(srange/slikezero_aud));        
    end

    like_vis = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_vis, srange), sigmasprime_vis).^2), sigmasprime_vis)/sqrt(2*pi);
    like_aud = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_aud, srange), sigmasprime_aud).^2), sigmasprime_aud)/sqrt(2*pi);
    like_vis = max(like_vis, NUMZERO);
    like_aud = max(like_aud, NUMZERO);
    
    % Compute prior, p(s)
    if priorinfo(3) == 1
        priorpdf = exp(-0.5*((srange - priorinfo(1))/priorinfo(2)).^2)/priorinfo(2)/sqrt(2*pi);
    else
        priorpdf = (priorinfo(3)*exp(-0.5*((priorinfo(1)-srange)/priorinfo(2)).^2)./priorinfo(2) + ...
            (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-srange)/priorinfo(5)).^2)./priorinfo(5))/sqrt(2*pi);
    end
    priorpdf = max(priorpdf, NUMZERO);
    priorpdf = priorpdf/(qtrapz(priorpdf, 1)*ds); % Normalize prior

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute optimal targets for vision and audition (C=2)
    theta_vis = [theta(1:6) theta(13:14)];
    sstaruni_vis = OptimalTargetUnimodal(X(X(:, 2)==1, 4), model, theta_vis, priorinfo, xrange_vis(:)', srange)';
    theta_aud = [theta(7:12) theta(13:14)];
    sstaruni_aud = OptimalTargetUnimodal(X(X(:, 2)==2, 3), model, theta_aud, priorinfo, xrange_aud(:)', srange);

    if nargout > 1 
        extras.sstaruni_vis = sstaruni_vis;    
        extras.sstaruni_aud = sstaruni_aud; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute strength of cue fusion

    % Compute marginal likelihood, p(x_vis, x_aud|C)

    % CASE C=2, Independent likelihoods
    likec2_vis = bsxfun(@times, priorpdf, like_vis);
    likec2_aud = bsxfun(@times, priorpdf, like_aud);
    likec2 = squeeze(bsxfun(@times, qtrapz(likec2_vis, 1)*ds, qtrapz(likec2_aud, 1)*ds)); 
        
    % CASE C=1, Likelihoods are not independent
    postpdfc1 = bsxfun(@times, priorpdf, bsxfun(@times, like_vis, like_aud));
    likec1 = qtrapz(postpdfc1, 1)*ds;    
    postpdfc1 = bsxfun(@rdivide, postpdfc1, likec1);
    
    %postpdfc1 = bsxfun(@times, priorpdf, bsxfun(@times, like_vis, like_aud));
    %likec1 = qtrapz(postpdfc1, 1)*ds;    
    %postpdfc1 = bsxfun(@rdivide, postpdfc1, likec1);

    
    % Compute weight for cue fusion
    switch model(15)
        case 1 % Model weight is Bayesian posterior p(C=1|x_1, x_2)
            likec1 = squeeze(likec1);
            postc1 = likec1*priorc1./(likec1*priorc1 + likec2*(1-priorc1));
            if priorc1cat ~= priorc1
                postc1cat = likec1*priorc1cat./(likec1*priorc1cat + likec2*(1-priorc1cat));
            end
            if softmaxbeta > 0
                postc1beta = 1./(1+exp(softmaxbeta*(1-2*postc1)));
            end
        case {2, 3} % Non-Bayesian, distance on x
            if softmaxbeta > 0
                postc1 = 1./(1 + exp(softmaxbeta_x*(abs(bsxfun(@minus, xrange_vis, xrange_aud)) - kcommon)));
                postc1beta = postc1;
            else
                postc1 = abs(bsxfun(@minus, xrange_vis, xrange_aud)) < kcommon;
            end
            postc1 = squeeze(postc1);
        case {4, 5} % Non-Bayesian, distance on shat
            if softmaxbeta > 0
                postc1 = 1./(1 + exp(softmaxbeta_x*(abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) - kcommon)));
            else
                postc1 = abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) < kcommon;                
            end
            % postc1 = squeeze(postc1);
        case 6 % Model weight is Bayesian posterior but extra stuff for unity
            likec1 = squeeze(likec1);
            postc1 = likec1*priorc1./(likec1*priorc1 + likec2*(1-priorc1));
            extrapostc1 = abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) < kcommon;
            extrapostc1 = squeeze(extrapostc1);
    end
    if nargout > 1 % Save variables for debug or data generation
        extras.postc1 = postc1;
        extras.postpdfc1 = postpdfc1;
        if priorc1cat ~= priorc1
            extras.postc1cat = postc1cat;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute optimal target for bimodal (C=1)

    if isinf(kappa)
        if sigmaloss == 0
            [~, index] = max(postpdfc1, [], 1);
            sstarbim(:, :) = srange(squeeze(index));        
        elseif isinf(sigmaloss)
            % sstarbim = squeeze(sum(bsxfun(@times, srange, postpdfc1), 1))*ds;
            temp = qtrapz(bsxfun(@times, srange, postpdfc1), 1);
            sstarbim(:, :) = temp(1, :, :)*ds;
        end

        if nargout > 1 % Save variables for debug or data generation
            extras.sstarbim = sstarbim;
        end
    end
        
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood for Gaussian likelihoods

    % Compute sensory likelihood std for vision
    if isinf(slikezero_vis) || slikezero_vis == 0
        sigmasprime2_vis = sigmalikezero_vis^2;
    elseif slikezero_vis > 0
        sigmasprime2_vis = sigmalikezero_vis^2*(1 + (xrange_vis/slikezero_vis).^2);    
    else % Negative szero for square-abs formula
        sigmasprime2_vis = sigmalikezero_vis^2*(1 + abs(xrange_vis/slikezero_vis)).^2;        
    end

    % Compute sensory likelihood std for audition
    if isinf(slikezero_aud) || slikezero_aud == 0
        sigmasprime2_aud = sigmalikezero_aud^2;
    elseif slikezero_aud > 0
        sigmasprime2_aud = sigmalikezero_aud^2*(1 + (xrange_aud/slikezero_aud).^2);    
    else % Negative szero for square-abs formula
        sigmasprime2_aud = sigmalikezero_aud^2*(1 + abs(xrange_aud/slikezero_aud)).^2;        
    end

    SSCALE = 8;
    srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
    ds = diff(srange(1:2));

    responsepdf = NaN(nTrials, 1);

    if nargout > 1 % Save variables for debug or data generation
        extras.xrange_vis = xrange_vis;
        extras.xrange_aud = xrange_aud;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute optimal targets for vision and audition (C=2)

    theta_vis = [theta(1:6) theta(13:14)];
    sstaruni_vis = OptimalTargetUnimodalGaussian(X(X(:, 2)==1, 4), model, theta_vis, priorinfo, xrange_vis(:)')';
    theta_aud = [theta(7:12) theta(13:14)];
    sstaruni_aud = OptimalTargetUnimodalGaussian(X(X(:, 2)==2, 3), model, theta_aud, priorinfo, xrange_aud(:)');
    
    if nargout > 1 
        extras.sstaruni_vis = sstaruni_vis;    
        extras.sstaruni_aud = sstaruni_aud; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute strength of cue fusion

    % Compute marginal likelihood, p(x_vis, x_aud|C)

    % CASE C=2, Independent likelihoods
    likec2_vis = priorinfo(3)*normpdf(xrange_vis, priorinfo(1), sqrt(priorinfo(2)^2 + sigmasprime2_vis)) + ...
        (1-priorinfo(3))*normpdf(xrange_vis, priorinfo(4), sqrt(priorinfo(5)^2 + sigmasprime2_vis));
    likec2_aud = priorinfo(3)*normpdf(xrange_aud, priorinfo(1), sqrt(priorinfo(2)^2 + sigmasprime2_aud)) + ...
        (1-priorinfo(3))*normpdf(xrange_aud, priorinfo(4), sqrt(priorinfo(5)^2 + sigmasprime2_aud));
    likec2 = squeeze(bsxfun(@times, likec2_vis, likec2_aud));

    % CASE C=1, Likelihoods are not independent
    sigma2sum = bsxfun(@plus, sigmasprime2_vis, sigmasprime2_aud);
    sigmatilde2 = bsxfun(@times, sigmasprime2_vis, sigmasprime2_aud)./sigma2sum;
    mutilde = bsxfun(@plus, bsxfun(@times, xrange_vis, sigmasprime2_aud), bsxfun(@times, xrange_aud, sigmasprime2_vis))./sigma2sum;
            
    %z1 = priorinfo(3) * bsxfun_normpdf(xrange_vis, xrange_aud, sqrt(sigma2sum)) .* ...
    %    bsxfun_normpdf(mutilde, priorinfo(1), sqrt(priorinfo(2)^2 + sigmatilde2));
    %z2 = (1-priorinfo(3)) * bsxfun_normpdf(xrange_vis, xrange_aud, sqrt(sigma2sum)) .* ...
    %    bsxfun_normpdf(mutilde, priorinfo(4), sqrt(priorinfo(5)^2 + sigmatilde2));
    
    logz1 = log(priorinfo(3)) + bsxfun_normlogpdf(xrange_vis, xrange_aud, sqrt(sigma2sum)) + ...
        bsxfun_normlogpdf(mutilde, priorinfo(1), sqrt(priorinfo(2)^2 + sigmatilde2));
    logz2 = log(1-priorinfo(3)) + bsxfun_normlogpdf(xrange_vis, xrange_aud, sqrt(sigma2sum)) + ...
        bsxfun_normlogpdf(mutilde, priorinfo(4), sqrt(priorinfo(5)^2 + sigmatilde2));
    z1 = exp(logz1); z2 = exp(logz2);
    
    % Compute weight for cue fusion
    switch model(15)
        case 1 % Model weight is Bayesian posterior p(C=1|x_1, x_2)
            likec1 = squeeze(z1 + z2);
            postc1 = likec1.*priorc1./(likec1.*priorc1 + likec2.*(1-priorc1));
            if priorc1cat ~= priorc1
                postc1cat = likec1*priorc1cat./(likec1*priorc1cat + likec2*(1-priorc1cat));
            end
            if softmaxbeta > 0
                postc1beta = 1./(1+exp(softmaxbeta*(1-2*postc1)));
            end            
        case {2, 3} % Non-Bayesian, distance on x    
            if softmaxbeta > 0
                postc1 = 1./(1 + exp(softmaxbeta_x*(abs(bsxfun(@minus, xrange_vis, xrange_aud)) - kcommon)));
                postc1beta = postc1;
            else
                postc1 = abs(bsxfun(@minus, xrange_vis, xrange_aud)) < kcommon;
            end
            postc1 = squeeze(postc1);
        case {4, 5} % Non-Bayesian, distance on shat
            postc1 = abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) < kcommon;
        case 6 % Model weight is Bayesian posterior but extra stuff for unity
            likec1 = squeeze(z1 + z2);
            postc1 = likec1.*priorc1./(likec1.*priorc1 + likec2.*(1-priorc1));
            extrapostc1 = abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) < kcommon;
            extrapostc1 = squeeze(extrapostc1);
    end
    if nargout > 1 % Save variables for debug or data generation
        extras.postc1 = postc1;
        if priorc1cat ~= priorc1
            extras.postc1cat = postc1cat;
        end    
    end
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute optimal target for bimodal (C=1)

    nz = max(logz1, logz2);
    logz1 = logz1 - nz;
    logz2 = logz2 - nz;    
    if isinf(kappa)
        if sigmaloss == 0
            error('MAP decision rule not supported for bimodal data yet.');
            %[~, index] = max(postpdfc1, [], 1);
            %sstarbim(:, :) = srange(squeeze(index));        
        elseif isinf(sigmaloss)
            mustar1 = (mutilde*priorinfo(2)^2 + priorinfo(1)*sigmatilde2)./(sigmatilde2 + priorinfo(2)^2);
            mustar2 = (mutilde*priorinfo(5)^2 + priorinfo(4)*sigmatilde2)./(sigmatilde2 + priorinfo(5)^2);
            sstarbim = squeeze((mustar1.*exp(logz1) + mustar2.*exp(logz2))./(exp(logz1) + exp(logz2)));
            % sstarbim = squeeze((mustar1.*z1 + mustar2.*z2)./(z1 + z2));
        end

        if nargout > 1 % Save variables for debug or data generation
            extras.sstarbim = sstarbim;
            extras.postpdfc1 = [];            
        end
    end    
    
end

% Clean up to prevent OUT OF MEMORY errors
clear postpdfc1 like_vis like_aud likec1 likec2 likec2_aud likec2_vis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute overall optimal targets

xpdf_nf_sum = zeros(size(X,1),1);
for iTask = 1:3
    f = (X(:, 2) == iTask);
    xpdf_vis = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_vis, alpha_vis*X(f, 4)), alpha_vis*sigmas_vis(f)).^2), alpha_vis*sigmas_vis(f))/sqrt(2*pi);
    nz_vis_trapz = qtrapz(xpdf_vis, 2);
    xpdf_vis = bsxfun(@rdivide, xpdf_vis, nz_vis_trapz*dx_vis); 
    xpdf_aud = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_aud, alpha_aud*X(f, 3)), alpha_aud*sigmas_aud(f)).^2), alpha_aud*sigmas_aud(f))/sqrt(2*pi);    
    nz_aud_trapz = qtrapz(xpdf_aud, 3);
    xpdf_aud = bsxfun(@rdivide, xpdf_aud, nz_aud_trapz*dx_aud);
    % Correction factor for SUM of xpdf (as opposed to QTRAPZ)
    xpdf_nf_sum(f) = nz_vis_trapz.*nz_aud_trapz./sum(xpdf_vis, 2)./sum(xpdf_aud,3);
    
    if model(14) == 3 || any(model(16) == [3 4]) % For probability matching use old method
        xpdf2{iTask} = bsxfun(@times, xpdf_vis, xpdf_aud); % These guys are already normalized
        % xpdf2{iTask} = bsxfun(@rdivide, xpdf2{iTask}, sum(sum(xpdf2{iTask}, 2),3)*dx_vis*dx_aud);
    else
        temp = bsxfun(@times, xpdf_vis, xpdf_aud);
        xpdf2{iTask} = reshape(temp, size(xpdf_vis, 1), XGRID*XGRID);        
    end
end

% Cycle across tasks
for iTask = 1:3
    f = (X(:, 2) == iTask);
    if sum(f) == 0; continue; end
        
    switch iTask
        case 1 % Visual localization
            % Deterministic decision making
            if isinf(kappa)                
                % Causal inference model for cue combination
                switch model(14)
                    case 1 % Model averaging
                        sstar(1, :, :) = bsxfun(@times, postc1, sstarbim) + bsxfun(@times, 1-postc1, sstaruni_vis);
                        responsepdf(f) = ComputeBimodalResponsePDF();
                        if nargout > 1; extras.sstar_bim_vis = squeeze(sstar(1,:,:)); end
                    case 2 % Model selection
                        if softmaxbeta > 0
                            sstar(1, :, :) = bsxfun(@times, postc1beta, sstarbim) + bsxfun(@times, 1-postc1beta, sstaruni_vis);
                        else
                            sstar(1, :, :) = bsxfun(@times, (postc1 >= rho), sstarbim) + bsxfun(@times, postc1 < rho, sstaruni_vis);
                        end
                        responsepdf(f) = ComputeBimodalResponsePDF();
                        if nargout > 1; extras.sstar_bim_vis = squeeze(sstar(1,:,:)); end
                    case 3 % Model probability matching                        
                        sstar1(1, :, :) = sstarbim;
                        sstar2(1, :, 1) = sstaruni_vis;
                        postc1shift(1, :, :) = postc1;
                        if wr == 0                            
                            responsepdf(f) = qtrapz(qtrapz(bsxfun(@times, xpdf2{iTask}, ...
                                bsxfun(@times, postc1shift, bsxfun_normpdf(sstar1, X(f, end), sigmar)) ...
                                + bsxfun(@times, 1 - postc1shift, bsxfun_normpdf(sstar2, X(f, end), sigmar))) ...
                                , 2), 3)*dx_vis*dx_aud;
                        else
                            sigmar1 = sqrt(sigmar^2 + wr^2*sstar1.^2);
                            sigmar2 = sqrt(sigmar^2 + wr^2*sstar2.^2);                            
                            responsepdf(f) = qtrapz(qtrapz(bsxfun(@times, xpdf2{iTask}, ...
                                bsxfun(@times, postc1shift, ...
                                bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus, sstar1, X(f, end)),sigmar1).^2),sigmar1) ) ...
                                + bsxfun(@times, 1 - postc1shift, ...
                                bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus, sstar2, X(f, end)),sigmar2).^2),sigmar2) ) )...
                                , 2), 3)*dx_vis*dx_aud/sqrt(2*pi);                            
                        end
                end
            end

        case 2 % Audio localization
            % Deterministic decision making
            if isinf(kappa)

                % Causal inference model for cue combination
                switch model(14)
                    case 1 % Model averaging
                        sstar(1, :, :) = bsxfun(@times, postc1, sstarbim) + bsxfun(@times, 1-postc1, sstaruni_aud);                        
                        responsepdf(f) = ComputeBimodalResponsePDF();
                        if nargout > 1; extras.sstar_bim_aud = squeeze(sstar(1,:,:)); end
                    case 2 % Model selection
                        if softmaxbeta > 0
                            sstar(1, :, :) = bsxfun(@times, postc1beta, sstarbim) + bsxfun(@times, 1-postc1beta, sstaruni_aud);
                        else
                            sstar(1, :, :) = bsxfun(@times, (postc1 >= rho), sstarbim) + bsxfun(@times, postc1 < rho, sstaruni_aud);
                        end
                        responsepdf(f) = ComputeBimodalResponsePDF();
                        if nargout > 1; extras.sstar_bim_aud = squeeze(sstar(1,:,:)); end
                    case 3 % Model probability matching
                        sstar1(1, :, :) = sstarbim;
                        sstar3(1, 1, :) = sstaruni_aud;
                        postc1shift(1, :, :) = postc1;
                        if wr == 0
                            responsepdf(f) = qtrapz(qtrapz(bsxfun(@times, xpdf2{iTask}, ...
                                bsxfun(@times, postc1shift, exp(-0.5*(bsxfun(@minus, sstar1, X(f, end))/sigmar).^2)) ...
                                + bsxfun(@times, 1 - postc1shift, exp(-0.5*(bsxfun(@minus, sstar3, X(f, end))/sigmar).^2))) ...
                                , 2), 3)*dx_vis*dx_aud/sqrt(2*pi)/sigmar;
                        else
                            sigmar1 = sqrt(sigmar^2 + wr^2*sstar1.^2);
                            sigmar3 = sqrt(sigmar^2 + wr^2*sstar3.^2);                            
                            responsepdf(f) = qtrapz(qtrapz(bsxfun(@times, xpdf2{iTask}, ...
                                bsxfun(@times, postc1shift, bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus, sstar1, X(f, end)),sigmar1).^2),sigmar1)) ...
                                + bsxfun(@times, 1 - postc1shift, bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus, sstar3, X(f, end)),sigmar3).^2),sigmar3)) ) ...
                                , 2), 3)*dx_vis*dx_aud/sqrt(2*pi);                            
                        end
                end
            end
            
        case 3 % Categorical response
            switch model(16)
                case {1, 2} % Standard threshold
                    if model(15) == 6
                        sstar(1, :, :) = 2 - ((postc1(:,:) + (1-postc1(:,:)).*extrapostc1(:,:)) > rho);                        
                    elseif priorc1cat ~= priorc1
                        sstar(1, :, :) = 2 - (postc1cat(:,:) > rho);                        
                    else
                        sstar(1, :, :) = 2 - (postc1(:,:) > rho);
                    end
                    sstarlin = logical(2 - sstar(:)');

                    f_one = X(:, end) == 1;
                    f_two = X(:, end) == 2;
                    
                    if 0
                        temp = bsxfun(@times,sum(xpdf2{iTask}(:, sstarlin), 2),xpdf_nf_sum(f))*dx_vis*dx_aud;

                        responsepdf(f & f_one) = temp(X(f, end) == 1);
                        responsepdf(f & f_two) = 1 - temp(X(f, end) == 2);
                        % responsepdf(f) = qtrapz(bsxfun(@times, xpdf2{iTask}, bsxfun(@eq, sstarlin, X(f, end))), 2)*dx_vis*dx_aud;

                        %responsepdf(f) = sum(sum(bsxfun(@times, xpdf2{iTask}, ...
                        %    bsxfun(@eq, sstar, X(f, end))), 2), 3)*dx_vis*dx_aud;
                    
                    else                                                
                        temp = bsxfun(@times,sum(bsxfun(@times,xpdf2{iTask},postc1beta(:)'), 2),xpdf_nf_sum(f))*dx_vis*dx_aud;
                        
                        responsepdf(f & f_one) = temp(X(f, end) == 1);
                        responsepdf(f & f_two) = 1 - temp(X(f, end) == 2);                        
                    end
                    
                    
                case {3, 4} % Probability matching
                    if priorc1cat ~= priorc1
                        postc1shift(1, :, :) = postc1cat;
                    else
                        postc1shift(1, :, :) = postc1;
                    end
                    responsepdf(f) = qtrapz(qtrapz(bsxfun(@times, xpdf2{iTask}, ...
                        bsxfun(@times, postc1shift, bsxfun(@eq, 1, X(f, end))) + bsxfun(@times, 1 - postc1shift, bsxfun(@eq, 2, X(f, end)))) ...
                        , 2), 3)*dx_vis*dx_aud;
            end
            
            % Fix probabilities            
            responsepdf(f) = min(max(responsepdf(f),0),1);
  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute prior-dependent lapse

if lambda > 0
    % Lapse over localization responses
    f = (X(:, 2) == 1) | (X(:, 2) == 2);
    if priorinfo(3) == 1
        responsepdf(f) = lambda*exp(-0.5*((priorinfo(1)-X(f, end))/priorinfo(2)).^2)./priorinfo(2)./sqrt(2*pi) ...
            + (1-lambda)*responsepdf(f);
    else
        responsepdf(f) = lambda*(priorinfo(3)*exp(-0.5*((priorinfo(1)-X(f, end))/priorinfo(2)).^2)./priorinfo(2) + ...
            (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-X(f, end))/priorinfo(5)).^2)./priorinfo(5)) ./sqrt(2*pi) ...
            + (1-lambda)*responsepdf(f);
    end
    % Lapse over categorical responses
    f = (X(:, 2) == 3);
    responsepdf(f) = lambdacat*(priorc1cat.*(X(f, end) == 1) + (1-priorc1cat).*(X(f, end) == 2)) ...
        + (1-lambdacat)*responsepdf(f);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize log likelihood

responsepdf = FIXEDLAPSEPDF + (1-FIXEDLAPSEPDF)*responsepdf;
if nargout > 1
    extras.responsepdf = responsepdf;
end

if sumover
    loglike = sum(log(responsepdf));
    varargout{1} = loglike;
else
    varargout{1} = responsepdf;
end    
if nargout > 1; varargout{2} = extras; end

return;

    %NORMPDF Normal probability density function (pdf).
    function y = normpdf(x,mu,sigma)
        y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    end

    % COMPUTEBIMODALRESPONSEPDF Compute response pdf for bimodal estimation.
    function temppdf = ComputeBimodalResponsePDF()
        sstarlin = sstar(:)'; % Convert from 2d matrix to 1d array
        temppdf = zeros(sum(f), 1);
        
        window = 6*sqrt(sigmar^2 + 4*wr^2*srange.^2);
        for ii = 1:length(srange)
            f2 = find(abs(X(f, end) - srange(ii)) < (ds/2)); % Identify all trials that fall in the same response bin
            if isempty(f2); continue; end % If the bin is empty, go to next bin
            index = find(abs(sstarlin - srange(ii)) < window(ii)); % Pick all x-positions whose sstar(x) falls within 6 SIGMAR of the response bin
            if isempty(index); continue; end % If there are no x-positions, go to next bin            
            sinv = sstarlin(index); % Pick all sstar(x)
            xx = X(f, end); % Pick all relevant responses
            nf = xpdf_nf_sum(f);
            if wr == 0
                temppdf(f2) = sum(bsxfun(@times, xpdf2{iTask}(f2, index), exp(-0.5*(bsxfun(@minus, sinv, xx(f2))/sigmar).^2)), 2).*nf(f2)*dx_vis*dx_aud/sqrt(2*pi)/sigmar;
            else
                sigmar1 = sqrt(sigmar^2 + wr^2*sinv.^2);
                temppdf(f2) = sum(bsxfun(@times, xpdf2{iTask}(f2, index), bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus, sinv, xx(f2)),sigmar1).^2),sigmar1)), 2).*nf(f2)*dx_vis*dx_aud/sqrt(2*pi);                
            end
        end        
    end

end

% Compute optimal target for unimodal condition
function sstar = OptimalTargetUnimodal(uniX, model, unitheta, priorinfo, xrange, srange)

    % This is almost zero
    NUMZERO = 1e-100;

    % Take model parameters
    sigmalikezero = unitheta(3);
    slikezero = unitheta(4);
    sigmaloss = unitheta(7);
    kappa = unitheta(8);

    ds = srange(2)-srange(1);

    % responsepdf = NaN(nTrials, 1);

    % Compute likelihood
    if slikezero >= 0
        sigmasprime = sigmalikezero*sqrt(1 + (srange/slikezero).^2);
    else % Negative slikezero for cosine noise formula
        % sigmasprime = sigmalikezero*(1 + abs(srange/slikezero));
        sigmasprime = sigmalikezero*sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90))/slikezero^2);
    end

    % Compute prior        
    if priorinfo(3) == 1
        priorpdf = exp(-0.5*((srange - priorinfo(1))/priorinfo(2)).^2)/priorinfo(2)/sqrt(2*pi);
    else
        priorpdf = (priorinfo(3)*exp(-0.5*((priorinfo(1)-srange)/priorinfo(2)).^2)./priorinfo(2) + ...
                (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-srange)/priorinfo(5)).^2)./priorinfo(5))/sqrt(2*pi);
    end
    if priorinfo(6) > 0 % Uniform prior
        priorpdf = priorpdf*(1-priorinfo(6)) + priorinfo(6)*(srange >= -priorinfo(7) & srange <= priorinfo(7))/(2*priorinfo(7));    
    end
    priorpdf = max(priorpdf, NUMZERO);    
    priorpdf = priorpdf/(qtrapz(priorpdf, 1)*ds); % Normalize

    % Compute unnormalized posterior for each measurement x in XRANGE
    postpdf = bsxfun(@times, priorpdf, ...
        bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2), sigmasprime));

    % Deterministic decision making
    if isinf(kappa)
        if sigmaloss == 0
            [~, index] = max(postpdf, [], 1);
            sstar = srange(index)';
        elseif isinf(sigmaloss)
            sstar = qtrapz(bsxfun(@times, srange, postpdf), 1)./qtrapz(postpdf, 1);
        end
    else
        error('Non-deterministic decision making not supported.');
    end
end

% Compute optimal target for unimodal condition
function sstar = OptimalTargetUnimodalGaussian(uniX, model, unitheta, priorinfo, xrange)

    % Take model parameters
    sigmalikezero = unitheta(3);
    slikezero = unitheta(4);
    sigmaloss = unitheta(7);
    kappa = unitheta(8);

    % Compute likelihood
    if isinf(slikezero) || slikezero == 0
        sigmasprime2 = sigmalikezero^2;    
    elseif slikezero > 0
        sigmasprime2 = sigmalikezero^2.*(1 + (xrange/slikezero).^2);
    else % Negative slikezero for cosine noise formula
        sigmasprime2 = sigmalikezero^2*(1 + 2*(90/pi)^2*(1 - cos(xrange*pi/90))/slikezero^2);
        % sigmasprime2 = sigmalikezero^2*(1 + abs(xrange/slikezero)).^2;
    end

    % Deterministic decision making
    if isinf(kappa)
        if sigmaloss == 0 % MAP decision rule
            error('MAP rule not supported for bimodal trials.');

        elseif isinf(sigmaloss) % MEAN decision rule
            mustar1 = (priorinfo(1)*sigmasprime2 + xrange*priorinfo(2)^2)./(priorinfo(2)^2 + sigmasprime2);
            mustar2 = (priorinfo(4)*sigmasprime2 + xrange*priorinfo(5)^2)./(priorinfo(5)^2 + sigmasprime2);
            logz1 = log(priorinfo(3)) + normlogpdf(xrange, priorinfo(1), sqrt(priorinfo(2)^2 + sigmasprime2));
            logz2 = log(1-priorinfo(3)) + normlogpdf(xrange, priorinfo(4), sqrt(priorinfo(5)^2 + sigmasprime2));      
            nz = max(logz1, logz2);
            z1 = exp(logz1 - nz);
            z2 = exp(logz2 - nz); 
            sstar = (mustar1.*z1 + mustar2.*z2)./(z1 + z2);
        end
    else
        % Compute decision distribution for each x (power of the posterior)
        error('Probabilistic decision rule not supported for bimodal trials.');        
    end
    
    
end