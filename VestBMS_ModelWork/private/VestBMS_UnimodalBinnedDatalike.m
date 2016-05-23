% VESTBMS_UNIMODALBINNEDDATALIKE Calculate log likelihood of unimodal 
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
function varargout = VestBMS_UnimodalBinnedDatalike(X,model,theta,priorinfo,bincenters,respbincenters,MAXRNG,XGRID,SSCALE,sumover,randomize)

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
MAXSD = 3;

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

% Compute sensory noise std per trial
if szero >= 0
    sigmas = sigmazero*sqrt(1 + (bincenters/szero).^2);
else % Negative szero for cosine noise formula
    sigmas = sigmazero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(bincenters*pi/90))/szero^2);
    % sigmas = sigmazero*(1 + abs(bincenters/szero));  % square-abs formula
end

% Measurements
xrange = alpha*linspace(max(min(bincenters-MAXSD*sigmas),-MAXRNG), min(max(bincenters+MAXSD*sigmas), MAXRNG), XGRID);
dx = xrange(2)-xrange(1);

% Add random jitter to XRANGE to estimate error in the computation of the 
% log likelihood due to the current grid size
if randomize; xrange = xrange + dx*(rand()-0.5); end

xpdf = bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, alpha*bincenters), alpha*sigmas).^2), alpha*sigmas)/sqrt(2*pi);
xpdf = bsxfun(@rdivide, xpdf, qtrapz(xpdf, 2)*dx);

if sigmaloss == 0; SSCALE = max(SSCALE, 24); end
srange = linspace(-MAXRNG, MAXRNG, MAXRNG*SSCALE*2 + 1)';

% Compute prior
if ~closedformflag
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type
if gaussianflag; likerange = xrange; else likerange = srange; end

if slikezero == 0
    sigmasprime = sigmalikezero;
elseif slikezero > 0
    sigmasprime = sigmalikezero*sqrt(1 + (likerange/slikezero).^2);
else % Negative slikezero for cosine noise formula
    sigmasprime = sigmalikezero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(likerange*pi/90))/slikezero^2);
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
            sstar = simpson2(bsxfun(@times, srange, postpdf), 1)./simpson2(postpdf,1);
            %sstar = trapz(bsxfun(@times, srange, postpdf), 1)./trapz(postpdf,1);
        elseif sigmaloss == -Inf                    % MEDIAN
            % Linear interpolation of median position from cdf
            cdf = bsxfun(@rdivide,cumtrapz(postpdf, 1),trapz(postpdf, 1));
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

% Bin optimal estimates to response bin centers
binbounds = [-MAXRNG; respbincenters(1:end-1) + 0.5*diff(respbincenters); MAXRNG];    
for iBin = 1:length(respbincenters)
    sstar(sstar > binbounds(iBin) &  sstar <= binbounds(iBin+1)) = respbincenters(iBin);
end

% Compute stimulus/response probability matrix (last dimension is response)
prmat = zeros(length(bincenters), length(respbincenters));
for iBin = 1:length(respbincenters)
    %prmat(:, iBin) = trapz(bsxfun(@times, xpdf, bsxfun(@eq, sstar, respbincenters(iBin))), 2)*dx;
    prmat(:, iBin) = simpson2(bsxfun(@times, xpdf, bsxfun(@eq, sstar, respbincenters(iBin))), 2)*dx;
end

if nargout > 1 % Save variables
    extras.xrange = xrange;
    extras.srange = srange;
    extras.sstar = sstar;
end

% Add prior-dependent lapse to response probability matrix
if lambda > 0    
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
        priorpdf = priorpdf./trapz(priorpdf,1);
        priorcdf = cumtrapz(priorpdf);
        cdfbin = [0, interp1(srange, priorcdf, lapsebinbounds(2:end-1))', 1];
        lapsepdf = diff(cdfbin);     
    end
    prmat = bsxfun(@plus, lambda*lapsepdf, (1-lambda)*prmat);
end
    
prmat = FIXEDLAPSEPDF + (1-FIXEDLAPSEPDF)*prmat;
prmat = prmat(:);

if nargout > 1
    extras.responsepdf = prmat;
end

if sumover
    loglike = sum(X(:).*log(prmat));
    varargout{1} = loglike;
else
    varargout{1} = prmat.^X(:);
end    
if nargout > 1; varargout{2} = extras; end

return;

    %LINSPACE Linearly spaced vector.
    function y = linspace(d1, d2, n)
        y = [d1 + ((0:n-2).*(d2-d1)/(n-1)), d2];
    end

    %NORMPDF Normal probability density function (pdf).
    function y = normpdf(x,mu,sigma)
        y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);
    end

    %NORMCDF Normal cumulative distribution function (cdf).
    function p = normcdf(x,mu,sigma)
        p = 0.5 * erfc(-(x-mu) ./ sigma ./ sqrt(2));
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