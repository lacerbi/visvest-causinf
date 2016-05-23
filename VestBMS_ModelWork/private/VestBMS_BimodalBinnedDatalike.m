% CUEBMS_BIMODALDATALIKE Calculate log likelihood of bimodal dataset X
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
function varargout = CueBMS_BimodalBinnedDatalike(X,model,theta,priorinfo,bincenters,maxranges,XGRID,SSCALE,sumover,randomize)

% Program constants
if nargin < 7 || isempty(XGRID) || isnan(XGRID); XGRID = 201; end
if nargin < 8 || isempty(SSCALE); SSCALE = 8; end
if nargin < 9 || isempty(sumover); sumover = 1; end
if nargin < 10 || isempty(randomize); randomize = 0; end

% Use closed form expression for Gaussians (ignores screen bounds!)
closedformflag = 0;

% When integrating a Gaussian, go up to this SDs away
MAXSD = 3;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1.5e-6;

% 1 visual localization, 2 audio or vestibular localization, 3 unity
NTASKS = 3; % Number of tasks 

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

%sigmar = theta(19); % Motor/placement noise
%wr = theta(22); % Weber's fraction for motor/placement noise
priorc1cat = theta(20); % P(C=1) for categorization trials
lambdacat = theta(21); % Categorisation trials lapse
if isnan(lambdacat); lambdacat = lambda; end

rho = 0.5; % Criterion for model selection
softmaxbeta = 1000; % Softmax parameter for model selection
softmaxbeta_multx = 0.1; % Multiplier for softmax beta when applied to x

% SD of prior over Delta for correlated prior (0 if uncorrelated)
priorsigmadelta = priorinfo(8);

% Model selection parameter (meaning depends on model)
priorc1 = priorinfo(end-1);
kcommon = priorinfo(end);

MAXRNG = maxranges(1);
if length(maxranges) > 1; MAXDELTA = maxranges(2); end

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

% Bin centers is a column vector
bincenters_vis = bincenters{1};
bincenters_aud = bincenters{2};
respbincenters = bincenters{3};

% Compute sensory noise std per trial for vision
if isinf(szero_vis) || szero_vis == 0
    sigmas_vis = sigmazero_vis*ones(length(bincenters_vis), 1);
elseif szero_vis > 0
    sigmas_vis = sigmazero_vis*sqrt(1 + (bincenters_vis/szero_vis).^2);    
else % Negative szero for cosine noise formula
    sigmas_vis = sigmazero_vis.*sqrt(1 + 2*(90/pi)^2*(1 - cos(bincenters_vis*pi/90))/szero_vis^2);    
    % sigmas_vis = sigmazero_vis*(1 + abs(bincenters_vis/szero_vis));        
end

% Compute sensory noise std per trial for audition
if isinf(szero_aud) || szero_aud == 0
    sigmas_aud = sigmazero_aud*ones(length(bincenters_aud), 1);
elseif szero_aud > 0
    sigmas_aud = sigmazero_aud*sqrt(1 + (bincenters_aud/szero_aud).^2);    
else % Negative szero for cosine noise formula (horrible way, I know)
    sigmas_aud = sigmazero_aud.*sqrt(1 + 2*(90/pi)^2*(1 - cos(bincenters_aud*pi/90))/szero_aud^2);
    % sigmas_aud = sigmazero_aud*(1 + abs(bincenters_aud/szero_aud));
end

% Measurements
xrange_vis = zeros(1, XGRID, 1);
xrange_aud = zeros(1, 1, XGRID);
xrange_vis(1, :, 1) = alpha_vis*linspace(max(min(bincenters_vis-MAXSD*sigmas_vis),-MAXRNG), min(max(bincenters_vis+MAXSD*sigmas_vis), MAXRNG), XGRID);
xrange_aud (1, 1, :) = alpha_aud*linspace(max(min(bincenters_aud-MAXSD*sigmas_aud),-MAXRNG), min(max(bincenters_aud+MAXSD*sigmas_aud), MAXRNG), XGRID);
dx_vis = xrange_vis(1, 2, 1) - xrange_vis(1, 1, 1);
dx_aud = xrange_aud(1, 1, 2) - xrange_aud(1, 1, 1);

% Add random jitter to XRANGE's to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize
    xrange_vis = xrange_vis + dx_vis*(rand()-0.5);
    xrange_aud = xrange_aud + dx_aud*(rand()-0.5);
end

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
ds = diff(srange(1:2));

if nargout > 1 % Save variables for debug or data generation
    extras.xrange_vis = xrange_vis;
    extras.xrange_aud = xrange_aud;
    extras.srange = srange;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type

if gaussianflag
    likerange_vis = xrange_vis; 
    likerange_aud = xrange_aud; 
else
    likerange_vis = srange; 
    likerange_aud = srange; 
end

% Compute sensory likelihood std for vision
if isinf(slikezero_vis) || slikezero_vis == 0
    sigmasprime_vis = sigmalikezero_vis;
elseif slikezero_vis > 0
    sigmasprime_vis = sigmalikezero_vis*sqrt(1 + (likerange_vis/slikezero_vis).^2);    
    % sigmasprime_vis = min(sigmasprime_vis, sigmalikezero_vis*sqrt(1 + (45/slikezero_vis).^2));    
else % Negative szero for cosine noise formula
    sigmasprime_vis = sigmalikezero_vis.*sqrt(1 + 2*(90/pi)^2*(1 - cos(likerange_vis*pi/90))/slikezero_vis^2);
    % sigmasprime_vis = sigmalikezero_vis*(1 + abs(likerange_vis/slikezero_vis));        
end

% Compute sensory likelihood std for audition
if isinf(slikezero_aud) || slikezero_aud == 0
    sigmasprime_aud = sigmalikezero_aud;
elseif slikezero_aud > 0
    sigmasprime_aud = sigmalikezero_aud*sqrt(1 + (likerange_aud/slikezero_aud).^2);    
else % Negative szero for cosine noise formula
    sigmasprime_aud = sigmalikezero_aud.*sqrt(1 + 2*(90/pi)^2*(1 - cos(likerange_aud*pi/90))/slikezero_aud^2);
    % sigmasprime_aud = sigmalikezero_aud*(1 + abs(likerange_aud/slikezero_aud));
end


if ~gaussianflag || ~closedformflag
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood for non-Gaussian likelihoods

    like_vis = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus,xrange_vis,srange), sigmasprime_vis).^2), sigmasprime_vis)/sqrt(2*pi);
    like_aud = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus,xrange_aud,srange), sigmasprime_aud).^2), sigmasprime_aud)/sqrt(2*pi);

    % Compute prior, p(s)
    if priorinfo(3) == 1
        priorpdf = exp(-0.5*((srange - priorinfo(1))/priorinfo(2)).^2)/priorinfo(2)/sqrt(2*pi);
    else
        priorpdf = (priorinfo(3)*exp(-0.5*((priorinfo(1)-srange)/priorinfo(2)).^2)./priorinfo(2) + ...
            (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-srange)/priorinfo(5)).^2)./priorinfo(5))/sqrt(2*pi);
    end
    priorpdf = priorpdf/(qtrapz(priorpdf, 1)*ds); % Normalize prior


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if priorsigmadelta == 0 % Uncorrelated prior on s
    
        % Compute optimal targets for vision and audition (C=2)
        theta_vis = [theta(1:6) theta(13:14)];
        sstaruni_vis = OptimalTargetUnimodal(model,theta_vis,priorinfo,xrange_vis(:)',srange)';

        theta_aud = [theta(7:12) theta(13:14)];
        sstaruni_aud = OptimalTargetUnimodal(model,theta_aud,priorinfo,xrange_aud(:)',srange);

        if nargout > 1 
            extras.sstaruni_vis = sstaruni_vis;    
            extras.sstaruni_aud = sstaruni_aud; 
        end
        
    else % Correlated prior on s and Delta
        
        % Introduce Delta
        deltarange(1,1,1,:) = linspace(-MAXDELTA, MAXDELTA, MAXDELTA*SSCALE*2 + 1)';
        ddelta = deltarange(2) - deltarange(1);
        
        % Joint prior on s and Delta
        prior2pdf = bsxfun(@times, ...
            normpdf(srange,priorinfo(1),priorinfo(2)), ...
            normpdf(deltarange,0,priorsigmadelta));
        prior2pdf = prior2pdf/(qtrapz(qtrapz(prior2pdf,1),4)*ds*ddelta);        
        
        % s_vis and s_aud in the s, Delta formulation
        stwo_vis = bsxfun(@minus,srange,0.5*deltarange);
        stwo_aud = bsxfun(@plus,srange,0.5*deltarange);
        
        % Recompute sensory likelihood std for vision
        if isinf(slikezero_vis) || slikezero_vis == 0
            sigmasprimetwo_vis = sigmalikezero_vis;
        elseif slikezero_vis > 0
            sigmasprimetwo_vis = sigmalikezero_vis*sqrt(1 + (stwo_vis/slikezero_vis).^2);    
        else % Negative szero for cosine noise formula
            sigmasprimetwo_vis = sigmalikezero_vis.*sqrt(1 + 2*(90/pi)^2*(1 - cos(stwo_vis*pi/90))/slikezero_vis^2);
        end

        % Recompute sensory likelihood std for audition
        if isinf(slikezero_aud) || slikezero_aud == 0
            sigmasprimetwo_aud = sigmalikezero_aud;
        elseif slikezero_aud > 0
            sigmasprimetwo_aud = sigmalikezero_aud*sqrt(1 + (stwo_aud/slikezero_aud).^2);    
        else % Negative szero for cosine noise formula
            sigmasprimetwo_aud = sigmalikezero_aud.*sqrt(1 + 2*(90/pi)^2*(1 - cos(stwo_aud*pi/90))/slikezero_aud^2);
        end
        
        % Sensory likelihood in s, Delta parametrization
        liketwo_vis = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus,xrange_vis,stwo_vis), sigmasprimetwo_vis).^2), sigmasprimetwo_vis)/sqrt(2*pi);
        liketwo_aud = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus,xrange_aud,stwo_aud), sigmasprimetwo_aud).^2), sigmasprimetwo_aud)/sqrt(2*pi);        
        liketwo = bsxfun(@times, liketwo_vis, liketwo_aud);
        
        % Posterior p(s,Delta|x_vis,x_aud)
        postpdftwo = bsxfun(@times,prior2pdf,liketwo);
        likec2 = qtrapz(qtrapz(postpdftwo,1),4); % Normalization factor
        
        % Compute MEAN of the posterior for vision and audition
        if isinf(kappa) && sigmaloss == Inf
            if ~isempty(X{1})       % Compute s*_vis(x_vis,x_aud)
                sstaruni_vis(:,:) = qtrapz(bsxfunandsum(@times,@times,stwo_vis.*prior2pdf,liketwo_vis,liketwo_aud,1,'qtrapz'),4)./likec2;
            end
            if ~isempty(X{2})       % Compute s*_aud(x_vis,x_aud)
                sstaruni_aud(:,:) = qtrapz(bsxfunandsum(@times,@times,stwo_aud.*prior2pdf,liketwo_vis,liketwo_aud,1,'qtrapz'),4)./likec2;
            end
        else
            error('Uncorrelated prior supports only deterministic MEAN decision rule.');
        end
        
        clear stwo_vis stwo_aud liketwo liketwo_vis liketwo_aud postpdftwo;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute strength of cue fusion

    % Compute marginal likelihood, p(x_vis, x_aud|C)

    % CASE C=2, Independent likelihoods
    if priorsigmadelta == 0 % Uncorrelated prior
        likec2_vis = bsxfun(@times, priorpdf, like_vis);
        likec2_aud = bsxfun(@times, priorpdf, like_aud);
        likec2 = squeeze(bsxfun(@times, qtrapz(likec2_vis, 1)*ds, qtrapz(likec2_aud, 1)*ds)); 
    else
        likec2 = squeeze(likec2)*ds*ddelta;
    end

    % CASE C=1, Likelihoods are not independent
    postpdfc1 = bsxfun(@times, priorpdf, bsxfun(@times, like_vis, like_aud));
    likec1 = qtrapz(postpdfc1, 1)*ds;
    postpdfc1 = bsxfun(@rdivide, postpdfc1, likec1);

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
                postc1 = 1./(1 + exp(softmaxbeta_multx*softmaxbeta*(abs(bsxfun(@minus, xrange_vis, xrange_aud)) - kcommon)));
            else
                postc1 = abs(bsxfun(@minus, xrange_vis, xrange_aud)) < kcommon;
            end
            postc1 = squeeze(postc1);
        case {4, 5} % Non-Bayesian, distance on shat
            if softmaxbeta > 0
                postc1 = 1./(1 + exp(softmaxbeta_multx*softmaxbeta*(abs(bsxfun(@minus, sstaruni_vis, sstaruni_aud)) - kcommon)));
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

    error('Closed form computation not supported (yet).');
    
    SSCALE = 8;
    srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
    ds = diff(srange(1:2));

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
    z1 = priorinfo(3) * exp(-0.5*bsxfun(@minus, xrange_vis, xrange_aud).^2./sigma2sum)./sqrt(sigma2sum) .* ...
        exp(-0.5*(mutilde - priorinfo(1)).^2./(priorinfo(2)^2 + sigmatilde2))./sqrt(priorinfo(2)^2 + sigmatilde2) / (2*pi);
    z2 = (1-priorinfo(3)) * exp(-0.5*bsxfun(@minus, xrange_vis, xrange_aud).^2./sigma2sum)./sqrt(sigma2sum) .* ...
        exp(-0.5*(mutilde - priorinfo(4)).^2./(priorinfo(5)^2 + sigmatilde2))./sqrt(priorinfo(5)^2 + sigmatilde2) / (2*pi);            

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
                postc1 = 1./(1 + exp(softmaxbeta_multx*softmaxbeta*(abs(bsxfun(@minus, xrange_vis, xrange_aud)) - kcommon)));
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

    if isinf(kappa)
        if sigmaloss == 0
            error('MAP decision rule not supported for bimodal data yet.');
            %[~, index] = max(postpdfc1, [], 1);
            %sstarbim(:, :) = srange(squeeze(index));        
        elseif isinf(sigmaloss)
            mustar1 = (mutilde*priorinfo(2)^2 + priorinfo(1)*sigmatilde2)./(sigmatilde2 + priorinfo(2)^2);
            mustar2 = (mutilde*priorinfo(5)^2 + priorinfo(4)*sigmatilde2)./(sigmatilde2 + priorinfo(5)^2);
            sstarbim = squeeze((mustar1.*z1 + mustar2.*z2)./(z1 + z2));        
        end

        if nargout > 1 % Save variables for debug or data generation
            extras.sstarbim = sstarbim;
            extras.postpdfc1 = [];            
        end
    end    
    
end

% Clean up to prevent OUT OF MEMORY errors
clear postpdfc1 like_vis like_aud likec1 likec2 likec2_aud likec2_vis temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute overall optimal targets

xpdf_vis = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_vis, alpha_vis*bincenters_vis), alpha_vis*sigmas_vis).^2), alpha_vis*sigmas_vis)/sqrt(2*pi);
nz_vis_trapz = qtrapz(xpdf_vis, 2);
xpdf_vis = bsxfun(@rdivide, xpdf_vis, nz_vis_trapz*dx_vis);
xpdf_aud = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange_aud, alpha_aud*bincenters_aud), alpha_aud*sigmas_aud).^2), alpha_aud*sigmas_aud)/sqrt(2*pi);    
nz_aud_trapz = qtrapz(xpdf_aud, 3);
xpdf_aud = bsxfun(@rdivide, xpdf_aud, nz_aud_trapz*dx_aud);
% Correction factor for SUM of xpdf (as opposed to QTRAPZ)
xpdf_nf_sum = nz_vis_trapz.*nz_aud_trapz./sum(xpdf_vis, 2)./sum(xpdf_aud,3);

if model(14) == 3 || any(model(16) == [3 4]) % For probability matching use old method
    xpdf2 = bsxfun(@times, xpdf_vis, xpdf_aud); % These guys are already normalized
    % xpdf2 = bsxfun(@rdivide, xpdf2, sum(sum(xpdf2, 2),3)*dx_vis*dx_aud);
else
    temp = bsxfun(@times, xpdf_vis, xpdf_aud);
    xpdf2 = reshape(temp, size(xpdf_vis, 1), XGRID*XGRID);        
end

% Probability response matrices (stimulus-response)
prmat = zeros(size(X{1},1) + size(X{2},1),length(respbincenters));
prmat_cat = zeros(size(X{3},1),2);

% Cycle across tasks
for iTask = 1:NTASKS
    if isempty(X{iTask}); continue; end
    
    switch iTask
        case 1 % Visual localization
            f = (1:size(X{1},1))';
            % Deterministic decision making
            if isinf(kappa)                
                % Causal inference model for cue combination
                switch model(14)
                    case 1 % Model averaging
                        sstar(1, :, :) = bsxfun(@times, postc1, sstarbim) + bsxfun(@times, 1-postc1, sstaruni_vis);
                        prmat(f,:) = ComputeBimodalResponsePDF(bincenters_vis);
                        if nargout > 1; extras.sstar_bim_vis = squeeze(sstar(1,:,:)); end
                    case 2 % Model selection
                        if softmaxbeta > 0
                            sstar(1, :, :) = bsxfun(@times, postc1beta, sstarbim) + bsxfun(@times, 1-postc1beta, sstaruni_vis);
                        else
                            sstar(1, :, :) = bsxfun(@times, (postc1 >= rho), sstarbim) + bsxfun(@times, postc1 < rho, sstaruni_vis);
                        end
                        prmat(f,:) = ComputeBimodalResponsePDF(bincenters_vis);
                        if nargout > 1; extras.sstar_bim_vis = squeeze(sstar(1,:,:)); end
                    case 3 % Model probability matching
                        sstar1(1, :, :) = BinResponses(sstarbim);
                        if isvector(sstaruni_vis)
                            sstar2(1, :, 1) = BinResponses(sstaruni_vis);
                        else
                            sstar2(1, :, :) = BinResponses(sstaruni_vis);                            
                        end
                        postc1shift(1, :, :) = postc1;

                        prmat_vis = zeros(length(bincenters_vis), length(respbincenters));
                        for rrBin = 1:length(respbincenters)
                            temp = bsxfun(@times, xpdf2, ...
                                bsxfun(@times, postc1shift, sstar1 == respbincenters(rrBin)) + ...
                                bsxfun(@times, 1-postc1shift, sstar2 == respbincenters(rrBin)));
                            prmat_vis(:, rrBin) = qtrapz(qtrapz(temp,2),3)*dx_vis*dx_aud;                            
                        end
                        prmat(f,:) = prmat_vis;
                end
            end

        case 2 % Audio localization
            f = size(X{1},1) + (1:size(X{2},1))';
            % Deterministic decision making
            if isinf(kappa)

                % Causal inference model for cue combination
                switch model(14)
                    case 1 % Model averaging
                        sstar(1, :, :) = bsxfun(@times, postc1, sstarbim) + bsxfun(@times, 1-postc1, sstaruni_aud);
                        prmat(f,:) = ComputeBimodalResponsePDF(bincenters_aud);
                        if nargout > 1; extras.sstar_bim_aud = squeeze(sstar(1,:,:)); end
                    case 2 % Model selection
                        if softmaxbeta > 0
                            sstar(1, :, :) = bsxfun(@times, postc1beta, sstarbim) + bsxfun(@times, 1-postc1beta, sstaruni_aud);
                        else
                            sstar(1, :, :) = bsxfun(@times, (postc1 >= rho), sstarbim) + bsxfun(@times, postc1 < rho, sstaruni_aud);
                        end
                        prmat(f,:) = ComputeBimodalResponsePDF(bincenters_aud);
                        if nargout > 1; extras.sstar_bim_aud = squeeze(sstar(1,:,:)); end
                    case 3 % Model probability matching
                        sstar1(1, :, :) = BinResponses(sstarbim);
                        if isvector(sstaruni_aud)
                            sstar3(1, 1, :) = BinResponses(sstaruni_aud);
                        else
                            sstar3(1, :, :) = BinResponses(sstaruni_aud);                            
                        end
                        postc1shift(1, :, :) = postc1;
                        
                        prmat_aud = zeros(length(bincenters_aud), length(respbincenters));
                        for rrBin = 1:length(respbincenters)
                            temp = bsxfun(@times, xpdf2, ...
                                bsxfun(@times, postc1shift, sstar1 == respbincenters(rrBin)) + ...
                                bsxfun(@times, 1-postc1shift, sstar3 == respbincenters(rrBin)));
                            prmat_aud(:, rrBin) = qtrapz(qtrapz(temp,2),3)*dx_vis*dx_aud;                            
                        end
                        prmat(f,:) = prmat_aud;
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
                    
                    temp = bsxfun(@times,sum(xpdf2(:, sstarlin), 2),xpdf_nf_sum)*dx_vis*dx_aud;                                        
                    prmat_cat(:,1) = temp;      % Prob of responding 1
                    prmat_cat(:,2) = 1 - temp;  % Prob of responding 2
                    
                case {3, 4} % Probability matching
                    if priorc1cat ~= priorc1
                        postc1shift(1, :, :) = postc1cat;
                    else
                        postc1shift(1, :, :) = postc1;
                    end
                    temp = qtrapz(qtrapz(bsxfun(@times,xpdf2,postc1shift),2),3)*dx_vis*dx_aud;
                    prmat_cat(:,1) = temp(:);      % Prob of responding 1
                    prmat_cat(:,2) = 1 - temp(:);  % Prob of responding 2                    
            end            
    end
end

% Fix probabilities
prmat = min(max(prmat,0),1);
prmat_cat = min(max(prmat_cat,0),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute prior-dependent lapse

if lambda > 0
    error('Lapse not supported yet.');
    
    % Compute lapse bin bounds
    % lapsebinbounds = binbounds;
    lapsebinbounds = [respbincenters(1)-0.5*diff(respbincenters(1:2)); ...
        respbincenters(1:end-1) + 0.5*diff(respbincenters); ...
        respbincenters(end)+0.5*diff(respbincenters(end-1:end))]';    
    if priorinfo(3) == 1 && priorinfo(6) == 0 % Use single Gaussian CDF
        cdfbin = normcdf(lapsebinbounds,priorinfo(1),priorinfo(2));
        lapsepdf = diff(cdfbin)/(cdfbin(end)-cdfbin(1));
    else % Compute CDF
        error('Need to fix this if want to use separate lapse bounds.');
        priorpdf = priorpdf./qtrapz(priorpdf,1);
        priorcdf = cumtrapz(priorpdf);
        cdfbin = [0, interp1(srange, priorcdf, lapsebinbounds(2:end-1))', 1];
        lapsepdf = diff(cdfbin);     
    end
    prmat = bsxfun(@plus, lambda*lapsepdf, (1-lambda)*prmat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finalize log likelihood

prmat = FIXEDLAPSEPDF + (1-FIXEDLAPSEPDF)*prmat;
prmat_cat = FIXEDLAPSEPDF + (1-FIXEDLAPSEPDF)*prmat_cat;
prmat = prmat(:); 
prmat_cat = prmat_cat(:);

if nargout > 1
    extras.responsepdf = [prmat; prmat_cat];
end

xx = [X{1}; X{2}]; xx_cat = X{3};

if sumover
    loglike = sum(xx(:).*log(prmat)) + sum(xx_cat(:).*log(prmat_cat));
    varargout{1} = loglike;
else
    varargout{1} = [prmat.^xx(:); prmat_cat.^xx_cat(:)];
end
if nargout > 1; varargout{2} = extras; end

return;

    %NORMPDF Normal probability density function (pdf).
    function y = normpdf(x,mu,sigma)
        y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    end

    % COMPUTEBIMODALRESPONSEPDF Compute response pdf for bimodal estimation.
    function prmat_mod = ComputeBimodalResponsePDF(bincenters_mod)
        sstarlin = sstar(:)'; % Convert from 2d matrix to 1d array        
        sstarlin = BinResponses(sstarlin);

        % Compute stimulus/response probability matrix (last dimension is response)
        prmat_mod = zeros(length(bincenters_mod), length(respbincenters));
        if length(respbincenters) == 2
            temp = sum(xpdf2(:, sstarlin == respbincenters(1)), 2).*xpdf_nf_sum*dx_vis*dx_aud;
            prmat_mod(:,1) = temp(:);
            prmat_mod(:,2) = 1-temp(:);            
        else
            for rBin = 1:length(respbincenters)
                prmat_mod(:,rBin) = sum(xpdf2(:, sstarlin == respbincenters(rBin)), 2).*xpdf_nf_sum*dx_vis*dx_aud;
            end
        end
    end

    %BINRESPONSES Bin optimal estimates to bin centers
    function sstar_binned = BinResponses(sstar_in)
        % Bin optimal estimates to bin centers
        sstar_binned = zeros(size(sstar_in));
        binbounds = [-MAXRNG; respbincenters(1:end-1) + 0.5*diff(respbincenters); MAXRNG];    
        for rBin = 1:length(respbincenters)
            sstar_binned(sstar_in > binbounds(rBin) &  sstar_in <= binbounds(rBin+1)) = respbincenters(rBin);
        end        
    end

end

%OPTIMALTARGETUNIMODAL Compute optimal target for unimodal condition 
% (optimal targets are NOT binned here)
function sstar = OptimalTargetUnimodal(model,unitheta,priorinfo,xrange,srange)

    closedformflag = 0;

    % Take model parameters
    sigmalikezero = unitheta(3);
    slikezero = unitheta(4);
    sigmaloss = unitheta(7);
    kappa = unitheta(8);

    % Is the internal noise Gaussian?
    if ( (slikezero == 0 || isinf(slikezero)) ...
            || ((model(4) == 6 || model(4) == 8)  && model(5) == 6)) ...
                                && priorinfo(6) == 0
        gaussianflag = 1;
    else
        gaussianflag = 0;
        if model(4) == 6 || model(5) == 6; error('Unsupported measurement-based likelihood and/or uniform prior.'); end
    end
        
    % Increase SSCALE if using MAP
    if isinf(kappa) && sigmaloss == 0 && SSCALE < 24; SSCALE = 24; end

    % ds = srange(2)-srange(1);
    
    % Compute prior
    if ~closedformflag
        if priorinfo(3) == 1
            priorpdf = exp(-0.5*((srange - priorinfo(1))/priorinfo(2)).^2)/priorinfo(2)/sqrt(2*pi);
        else
            priorpdf = (priorinfo(3)*exp(-0.5*((priorinfo(1)-srange)/priorinfo(2)).^2)./priorinfo(2) + ...
                    (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-srange)/priorinfo(5)).^2)./priorinfo(5))/sqrt(2*pi);
        end
        if priorinfo(6) > 0 % Uniform prior
            priorpdf = priorpdf*(1-priorinfo(6)) + priorinfo(6)*(srange >= -priorinfo(7) & srange <= priorinfo(7))/(2*priorinfo(7));    
        end
        % priorpdf = priorpdf/(sum(priorpdf, 1)*ds); % Normalize
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute likelihood; variable range depends on type
    if gaussianflag; likerange = xrange; else likerange = srange; end

    if slikezero == 0
        sigmasprime = sigmalikezero;
    elseif slikezero > 0
        sigmasprime = sigmalikezero*sqrt(1 + (likerange/slikezero).^2);
    else % Negative slikezero for square-abs formula
        sigmasprime = sigmalikezero*(1 + abs(likerange/slikezero));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute optimal estimate

    if closedformflag && isinf(sigmaloss) %  Closed form (ignores boundaries)
        sigmasprime2 = sigmasprime.^2;
        mustar1 = (priorinfo(1)*sigmasprime2 + xrange*priorinfo(2)^2)./(priorinfo(2)^2 + sigmasprime2);
        mustar2 = (priorinfo(4)*sigmasprime2 + xrange*priorinfo(5)^2)./(priorinfo(5)^2 + sigmasprime2);
        z1 = priorinfo(3)*normpdf(xrange, priorinfo(1), sqrt(priorinfo(2)^2 + sigmasprime2));
        z2 = (1-priorinfo(3))*normpdf(xrange, priorinfo(4), sqrt(priorinfo(5)^2 + sigmasprime2));        
        sstar = (mustar1.*z1 + mustar2.*z2)./(z1 + z2);
    else
        % Compute unnormalized posterior for each measurement x in XRANGE
        postpdf = bsxfun(@times, priorpdf, ...
            bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2), sigmasprime));

        % Deterministic decision making
        if isinf(kappa)
            if sigmaloss == 0                           % MAP
                [~, index] = max(postpdf, [], 1);
                sstar = srange(index)';
            elseif sigmaloss == Inf                     % MEAN
                sstar = qtrapz(bsxfun(@times, srange, postpdf), 1)./qtrapz(postpdf, 1);
            elseif sigmaloss == -Inf                    % MEDIAN
                % Linear interpolation of median position from cdf
                cdf = bsxfun(@rdivide,cumtrapz(postpdf, 1),qtrapz(postpdf, 1));
                [~,pos] = max(cdf >= 0.5, [], 1);
                p0 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + max(pos - 1, 1));
                p1 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + pos);
                coeff1 = (1./(1 + (0.5 - p0)./(p1 - 0.5)))';
                sstar = ((1-coeff1).*srange(pos) + (coeff1).*srange(max(pos-1, 1)))';
            end
        else
            error('Stochastic decision making not supported.');
        end
    end
end