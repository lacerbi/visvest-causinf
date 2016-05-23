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
function varargout = CueBMS_UnimodalDatalike_fixedx(X,model,theta,priorinfo,MAXRNG,xrange,SSCALE,sumover,randomize)

% Program constants
if nargin < 7 || isempty(SSCALE); SSCALE = 8; end
if nargin < 8 || isempty(sumover); sumover = 1; end
if nargin < 9 || isempty(randomize); randomize = 0; end

persistent infostruct;

if isempty(infostruct)
    infostruct.xrange = randn(1000,1);
end

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1.5e-6;

% This is almost zero
NUMZERO = 1e-80;

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

% Increase SSCALE if using MAP
if isinf(kappa) && sigmaloss == 0 && SSCALE < 24; SSCALE = 24; end

% Compute sensory noise std per trial
if szero >= 0
    sigmas = sigmazero*sqrt(1 + (X(:, 2)/szero).^2);
else % Negative szero for cosine noise formula
    % sigmas = sigmazero*(1 + abs(X(:,2)/szero));
    sigmas = sigmazero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(X(:,2)*pi/90))/szero^2);    
end

% Measurements
xrange = alpha*infostruct.xrange(1:nTrials)';

srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';
ds = srange(2)-srange(1);

% Compute likelihood
if slikezero >= 0
    sigmasprime = sigmalikezero*sqrt(1 + (srange/slikezero).^2);
    sigmatilde2_fun = @(s) sigmalikezero^2.*(1 + (s/slikezero).^2);
else % Negative slikezero for cosine noise formula
    % sigmasprime = sigmalikezero*(1 + abs(srange/slikezero));
    sigmasprime = sigmalikezero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90))/slikezero^2);
    sigmatilde2_fun = @(s) sigmalikezero^2.*(1 + 2*(90/pi)^2*(1 - cos(s*pi/90))/slikezero^2);
end

prior.w = [priorinfo(3),1-priorinfo(3)];
prior.mu = [priorinfo(1),priorinfo(4)];
prior.sigma = [priorinfo(2),priorinfo(5)];


% Compute prior
%priorpdf = computePrior(srange);
%priorpdf = priorpdf/(simpson2(priorpdf, 1)*ds); % Normalize

% Compute unnormalized posterior for each measurement x in XRANGE
%postpdf = bsxfun(@times, priorpdf, ...
%    bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2), sigmasprime));



% Deterministic decision making
if isinf(kappa)
    if sigmaloss == 0                           % MAP
        [~, index] = max(postpdf, [], 1);
        sstar = srange(index)';
    elseif sigmaloss == Inf                     % MEAN
        %sstar = (simpson2(bsxfun(@times, srange, postpdf), 1)./simpson2(postpdf, 1))';
        [xrange_ord,ord] = sort(xrange);
        xpivot = xrange_ord([1:20:end-1,end]);
        sstar(ord,1) = shatest1d(xpivot,xrange_ord,sigmatilde2_fun,prior);        
        
    elseif sigmaloss == -Inf                    % MEDIAN
        % Linear interpolation of median position from cdf
        cdf = bsxfun(@rdivide,cumtrapz(postpdf, 1),trapz(postpdf, 1));
        [~,pos] = max(cdf >= 0.5, [], 1);
        p0 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + max(pos - 1, 1));
        p1 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + pos);
        coeff1 = (1./(1 + (0.5 - p0)./(p1 - 0.5)))';
        sstar = ((1-coeff1).*srange(pos) + (coeff1).*srange(max(pos-1, 1)))';
    end
    
    if nargout > 1
        extras.xrange = xrange;
        extras.srange = srange;
        extras.sstar = sstar;
    end
        
    if wr == 0; sigmartot = sigmar; else sigmartot = sqrt(sigmar^2 + wr^2*sstar.^2); end    
    xpdf = exp(-0.5*((xrange' - alpha*X(:, 2))./(alpha*sigmas)).^2)./(alpha*sigmas)/sqrt(2*pi);
    responsepdf = exp(-0.5*(bsxfun(@minus, sstar, X(:, 3))/sigmartot).^2)/sigmartot/sqrt(2*pi).*xpdf;
else
    error('Stochastic decision making not supported.'); 
end

lambda = 0;

% Add prior-dependent lapse
if lambda > 0
    lapsepdf = computePrior(X(:,3));
    responsepdf = lambda*lapsepdf + (1-lambda)*responsepdf;
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

    %COMPUTEPRIOR Compute unnormalized prior
    function y = computePrior(S)
        if priorinfo(3) == 1
            y = exp(-0.5*((S - priorinfo(1))/priorinfo(2)).^2)/priorinfo(2)/sqrt(2*pi);
        else
            y = (priorinfo(3)*exp(-0.5*((priorinfo(1)-S)/priorinfo(2)).^2)./priorinfo(2) + ...
                    (1 - priorinfo(3))*exp(-0.5*((priorinfo(4)-S)/priorinfo(5)).^2)./priorinfo(5))/sqrt(2*pi);
        end
        if priorinfo(6) > 0 % Uniform prior
            if isinf(priorinfo(6)) % Multiplied uniform
                y(S < -priorinfo(7) | S > priorinfo(7)) = eps;
            else % Mixture with uniform
                y = y*(1-priorinfo(6)) + priorinfo(6)*(S >= -priorinfo(7) & S <= priorinfo(7))/(2*priorinfo(7));
            end
        end
        y = max(y, NUMZERO);
    end

    %LINSPACE Linearly spaced vector.
    function y = linspace(d1, d2, n)
        y = [d1 + ((0:n-2).*(d2-d1)/(n-1)), d2];
    end

    %SIMPSON2 Integration with Simpson's rule along one direction on a 2D
    %matrix.
    function y = simpson2(x,dim)
        w = [17 59 43 49]/48;
        if dim == 1
            y = w(1)*x(1,:) + w(2)*x(2,:) + w(3)*x(3,:) + w(4)*x(4,:) + ...
                sum(x(5:end-4,:),1) + ...
                w(4)*x(end-3,:) + w(3)*x(end-2,:) + w(2)*x(end-1,:) + w(1)*x(end,:);            
        else
            y = w(1)*x(:,1) + w(2)*x(:,2) + w(3)*x(:,3) + w(4)*x(:,4) + ...
                sum(x(:,5:end-4),2) + ...
                w(4)*x(:,end-3) + w(3)*x(:,end-2) + w(2)*x(:,end-1) + w(1)*x(:,end);            
        end
    end

end

