% CUEBMS_UNIMODALDATALIKE Calculate log likelihood of unimodal dataset X
% under model MODEL and parameter set THETA.
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
function varargout = CueBMS_UnimodalGaussianDatalike(X,model,theta,priorinfo,XGRID,sumover,randomize)

% Program constants
if nargin < 5 || isempty(XGRID) || isnan(XGRID); XGRID = 201; end
if nargin < 6 || isempty(sumover); sumover = 1; end
if nargin < 7 || isempty(randomize); randomize = 0; end

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;

% Screen bounds (in degrees)
MAXRNG = 45;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1.5e-6;

% Take model parameters
sigmazero = theta(1);
szero = theta(2);
sigmalikezero = theta(3);
slikezero = theta(4);
alpha = theta(5);
beta = theta(6);
sigmaloss = theta(7);
kappa = theta(8);
lambda = theta(9);
sigmar = theta(13); % Motor/placement noise
wr = theta(14); % Weber's fraction for motor/placement noise

% Number of trials in this condition
nTrials = size(X, 1);

% Compute sensory noise std per trial
if isinf(szero) || szero == 0
    sigmas = sigmazero;
elseif szero > 0
    sigmas = sigmazero*sqrt(1 + (X(:, 2)/szero).^2);    
else % Negative szero for cosine noise formula
    % sigmas = sigmazero*(1 + abs(X(:, 2)/szero));
    sigmas = sigmazero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(X(:,2)*pi/90))/szero^2);    
end

% Measurements
xrange = alpha*linspace(max(min(X(:,2)-MAXSD*sigmas),-MAXRNG), min(max(X(:,2)+MAXSD*sigmas,MAXRNG)), XGRID);
dx = xrange(2)-xrange(1);

% Add random jitter to XRANGE to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize; xrange = xrange + dx*(rand()-0.5); end

xpdf = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, alpha*X(:, 2)), alpha*sigmas).^2), alpha*sigmas)/sqrt(2*pi);
xpdf = bsxfun(@rdivide, xpdf, qtrapz(xpdf, 2)*dx);

% Compute likelihood
if isinf(slikezero) || slikezero == 0
    sigmasprime2 = sigmalikezero^2;    
elseif slikezero > 0
    sigmasprime2 = sigmalikezero^2.*(1 + (xrange/slikezero).^2);
else % Negative slikezero for cosine noise formula
    % sigmasprime2 = sigmalikezero^2*(1 + abs(xrange/slikezero)).^2;
    sigmasprime2 = sigmalikezero^2.*(1 + 2*(90/pi)^2*(1 - cos(likerange*pi/90))/slikezero^2);
end

responsepdf = NaN(nTrials, 1);

% Deterministic decision making
if isinf(kappa)
    if sigmaloss == 0 % MAP decision rule
        SSCALE = 24;
        srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
        
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
        
        % Likelihood
        sigmasprime = sqrt(sigmasprime2);
        
        % Compute unnormalized posterior for each measurement x in XRANGE
        postpdf = bsxfun(@times, priorpdf, ...
            bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2), sigmasprime));
        
        % Compute MAP
        [~, index] = max(postpdf, [], 1);
        sstar = srange(index)';
        
    elseif isinf(sigmaloss) % MEAN decision rule
        if priorinfo(3) == 1
            sstar = (priorinfo(1)*sigmasprime2 + xrange*priorinfo(2)^2)./(priorinfo(2)^2 + sigmasprime2);
        else
            mustar1 = (priorinfo(1)*sigmasprime2 + xrange*priorinfo(2)^2)./(priorinfo(2)^2 + sigmasprime2);
            mustar2 = (priorinfo(4)*sigmasprime2 + xrange*priorinfo(5)^2)./(priorinfo(5)^2 + sigmasprime2);
            z2overz1 = (1-priorinfo(3))/priorinfo(3) * ...
                exp(-0.5*((xrange-priorinfo(4)).^2/(priorinfo(5)^2 + sigmasprime2) ...
                - (xrange-priorinfo(1)).^2/(priorinfo(2)^2 + sigmasprime2))) ...
                ./ sqrt(priorinfo(5)^2 + sigmasprime2) .* sqrt(priorinfo(2)^2 + sigmasprime2) ;
            z2overz1 = min(z2overz1, 1e80);
            sstar = (mustar1 + mustar2.*z2overz1)./(1 + z2overz1);
            %z1 = priorinfo(3)*normpdf(xrange, priorinfo(1), sqrt(priorinfo(2)^2 + sigmasprime2));
            %z2 = (1-priorinfo(3))*normpdf(xrange, priorinfo(4), sqrt(priorinfo(5)^2 + sigmasprime2));        
            %sstar = (mustar1.*z1 + mustar2.*z2)./(z1 + z2);
        end
    end
    
    if any(isnan(sstar) | isinf(sstar) | ~isreal(sstar))
        priorinfo(3)
        z2overz1
        sstar
    end
    
    
    if nargout > 1
        extras.xrange = xrange;
        extras.sstar = sstar;
    end
    
    if wr == 0
        responsepdf = qtrapz(bsxfun(@times, xpdf, ...
            exp(-0.5*(bsxfun(@minus, sstar, X(:, 3))/sigmar).^2)), 2)*dx/sigmar/sqrt(2*pi);
    else
        sigmartot = sqrt(sigmar^2 + wr^2*sstar.^2);
        responsepdf = qtrapz(bsxfun(@times, xpdf, ...
            bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, sstar, X(:, 3)), sigmartot).^2), sigmartot)), 2)*dx/sqrt(2*pi);
    end
else     
    % Compute decision distribution for each x (power of the posterior)
        
end

% Add prior-dependent lapse
if lambda > 0
    if priorinfo(3) == 1
        responsepdf = lambda*exp(-0.5*((priorinfo(1)-X(:, 3))/priorinfo(2)).^2)./priorinfo(2)./sqrt(2*pi) ...
            + (1-lambda)*responsepdf;
    else
        responsepdf = lambda*(priorinfo(3)*exp(-0.5*((priorinfo(1)-X(:, 3))/priorinfo(2)).^2)./priorinfo(2) + ...
            (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-X(:, 3))/priorinfo(5)).^2)./priorinfo(5)) ./sqrt(2*pi) ...
            + (1-lambda)*responsepdf;        
    end
end

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

end