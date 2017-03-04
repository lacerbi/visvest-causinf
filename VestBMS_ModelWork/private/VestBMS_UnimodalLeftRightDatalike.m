% VESTBMS_UNIMODALLEFTRIGHTDATALIKE Calculate (log) likelihood of unimodal 
% dataset X under model MODEL and parameter set THETA.
%
% X is a unimodal data matrix. Each row is a trial. For a given row, the 
% columns contain data for:
% X(1) Number of trial (unused),
% X(2) Stimulus position (deg),
% X(3) Response (deg).
%
% MODEL is the model class parameter vector.
% MODEL(1) is the number of effective parameters.
%
% THETA is the model-dependent parameter vector for the current condition. 
% The values of THETA are:
% THETA(1) sigmazero; THETA(2) w; 
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
function varargout = VestBMS_UnimodalLeftRightDatalike(X,model,theta,priorinfo,bincenters,MAXRNG,XGRID,SSCALE,sumover,randomize)

% Program constants
if nargin < 7 || isempty(XGRID) || isnan(XGRID); XGRID = 201; end
if nargin < 8 || isempty(SSCALE); SSCALE = 8; end
if nargin < 9 || isempty(sumover); sumover = 1; end
if nargin < 10 || isempty(randomize); randomize = 0; end

%SSCALE = 8;
%XGRID = 251;

% Use closed form expression for Gaussians (ignores screen bounds!)
closedformflag = 0;

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1e-4;

MAXRNG_XMEAS = 180;

% Take model parameters
sigmazero = theta(1);
w = theta(2);
sigmalikezero = theta(3);
wlike = theta(4);
alpha_rescaling = 1;
beta_softmax = theta(8);
lambda = theta(9);

% Is the internal noise Gaussian?
if ( wlike == 0 ...
        || ((model(4) == 6 || model(4) == 8)  && model(5) == 6))
    gaussianflag = 1;
else
    gaussianflag = 0;
    if model(4) == 6 || model(5) == 6; error('Unsupported measurement-based likelihood and/or uniform prior.'); end
end

% Skip the decision process if it is irrelevant
skipdecision = (priorinfo(1) == 0) && (beta_softmax == Inf);

% Compute sensory noise std per trial
if w >= 0; noisemodel = 'A'; else noisemodel = 'C'; end
sigmas = VestBMS_sensoryNoise(noisemodel,bincenters,sigmazero,w);
if isscalar(sigmas); sigmas = sigmas*ones(size(bincenters)); end

% Measurements
xrange = alpha_rescaling*linspace(max(min(bincenters-MAXSD*sigmas),-MAXRNG_XMEAS), min(max(bincenters+MAXSD*sigmas), MAXRNG_XMEAS), XGRID);
dx = xrange(2)-xrange(1);

% Wrap large noisy measurement around circle
if MAXRNG_XMEAS >= 180 && ...
        (min(bincenters-MAXSD*sigmas) <= -180 || max(bincenters+MAXSD*sigmas) >= 180)
    wraparound = 1;
else
    wraparound = 0;
end

% Add random jitter to XRANGE to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize; xrange = xrange + dx*(rand()-0.5); end

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';

if skipdecision    
    prright = bsxfun_normcdf(alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters, alpha_rescaling*sigmas) - ...
        bsxfun_normcdf(0, alpha_rescaling*bincenters, alpha_rescaling*sigmas);
    if wraparound
        prright = prright + bsxfun_normcdf(alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters + 360, alpha_rescaling*sigmas) - ...
            bsxfun_normcdf(0, alpha_rescaling*bincenters + 360, alpha_rescaling*sigmas);
        prright = prright + bsxfun_normcdf(alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters - 360, alpha_rescaling*sigmas) - ...
            bsxfun_normcdf(0, alpha_rescaling*bincenters - 360, alpha_rescaling*sigmas);
    end
    
    prleft = bsxfun_normcdf(0, alpha_rescaling*bincenters, alpha_rescaling*sigmas) - ...
        bsxfun_normcdf(-alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters, alpha_rescaling*sigmas);
    if wraparound
        prleft = prleft + bsxfun_normcdf(0, alpha_rescaling*bincenters + 360, alpha_rescaling*sigmas) - ...
            bsxfun_normcdf(-alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters + 360, alpha_rescaling*sigmas);
        prleft = prleft + bsxfun_normcdf(0, alpha_rescaling*bincenters - 360, alpha_rescaling*sigmas) - ...
            bsxfun_normcdf(-alpha_rescaling*MAXRNG_XMEAS, alpha_rescaling*bincenters - 360, alpha_rescaling*sigmas);
    end
    
    % Compute stimulus/response probability matrix (last dimension is response)
    prmat = zeros(numel(bincenters), 2);
    prmat(:, 1) = prleft./(prright + prleft);
    prmat(:, 2) = prright./(prright + prleft);
    
else

    % Compute prior
    if ~closedformflag
        priorpdf = bsxfun_normpdf(srange,priorinfo(1),priorinfo(2));
    end

    % Compute likelihood; variable range depends on type
    if gaussianflag; likerange = xrange; else likerange = srange; end
    if wlike >= 0; likemodel = 'A'; else likemodel = 'C'; end
    sigmasprime = VestBMS_sensoryNoise(likemodel,likerange,sigmalikezero,wlike);

    if wraparound
        like_meas = bsxfun_normpdf(xrange, srange, sigmasprime) + ...
            bsxfun_normpdf(xrange, srange + 360, sigmasprime) + ...
            bsxfun_normpdf(xrange, srange - 360, sigmasprime);
            % The fact that the first and last column of XRANGE is counted
            % twice, due to wrapping, is accounted for by the subsequent QTRAPZ
    else
        like_meas = bsxfun_normpdf(xrange, srange, sigmasprime);
    end

    % Compute unnormalized posterior for each measurement x in XRANGE
    postpdf = bsxfun(@times, priorpdf, like_meas);
    
    if isinf(beta_softmax) && 0
        % Deterministic decision making
        error('Deterministic decision making not supported.');

        % Linear interpolation of median position from cdf
        cdf = bsxfun(@rdivide,cumtrapz(postpdf, 1),trapz(postpdf, 1));
        [~,pos] = max(cdf >= 0.5, [], 1);
        p0 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + max(pos - 1, 1));
        p1 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + pos);
        coeff1 = (1./(1 + (0.5 - p0)./(p1 - 0.5)))';
        sstar = ((1-coeff1).*srange(pos) + (coeff1).*srange(max(pos-1, 1)))';

        % Bin optimal estimates to response bin centers
        %binbounds = [-MAXRNG; respbincenters(1:end-1) + 0.5*diff(respbincenters); MAXRNG];    
        %for iBin = 1:length(respbincenters)
        %    sstar(sstar > binbounds(iBin) &  sstar <= binbounds(iBin+1)) = respbincenters(iBin);
        %end

    else
        % Compute observer's posterior probability of rightward motion
        postright = VestBMS_PostRight(postpdf);

        % Probability of rightward response
        prright = 1./(1 + ((1-postright)./postright).^beta_softmax);
    end
    
    % Compute noise distribution
    if wraparound
        xpdf = bsxfun_normpdf(xrange, alpha_rescaling*bincenters, alpha_rescaling*sigmas) + ...
            bsxfun_normpdf(xrange, alpha_rescaling*bincenters - 360, alpha_rescaling*sigmas) + ...
            bsxfun_normpdf(xrange, alpha_rescaling*bincenters + 360, alpha_rescaling*sigmas);
    else
        xpdf = bsxfun_normpdf(xrange, alpha_rescaling*bincenters, alpha_rescaling*sigmas);
    end
    xpdf = bsxfun(@rdivide, xpdf, qtrapz(xpdf, 2)*dx);

    % Compute stimulus/response probability matrix (last dimension is response)
    prmat = zeros(numel(bincenters), 2);
    prmat(:, 2) = simpson1(bsxfun(@times, xpdf, prright), 2)*dx;
    prmat(:, 1) = 1 - prmat(:, 2);
    
end

% Finalize log likelihood
prmat = lambda/2 + (1-lambda)*prmat;
prmat = 0.5*FIXEDLAPSEPDF + (1-FIXEDLAPSEPDF)*prmat;
prmat = prmat(:);

if nargout > 1
    extras.xrange = xrange;
    extras.srange = srange;
    extras.responsepdf = prmat;
end

if sumover(1)
    loglike = sum(X(:).*log(prmat));
    varargout{1} = loglike;
else    
    varargout{1} = loglikmat2vec(log(prmat),X(:));
    % varargout{1} = prmat.^X(:);
end    
if nargout > 1; varargout{2} = extras; end

end