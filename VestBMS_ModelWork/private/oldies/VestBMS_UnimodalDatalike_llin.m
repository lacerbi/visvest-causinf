%CUEBMS_UNIMODALDATALIKE_LLIN Log likelihood of uni-sensory dataset with 
%  piecewise linear approximation.
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
function varargout = CueBMS_UnimodalDatalike_llin(X,model,theta,priorinfo,MAXRNG,XGRID,SSCALE,sumover,randomize)

% Program constants
if nargin < 6 || isempty(XGRID) || isnan(XGRID); XGRID = 201; end
if nargin < 7 || isempty(SSCALE); SSCALE = 8; end
if nargin < 8 || isempty(sumover); sumover = 1; end
if nargin < 9 || isempty(randomize); randomize = 0; end

SSCALE = 32;

% Use taylor approximation for the optimal estimate
dotaylor = 0;

% Use piecewise linear approximation for the probability of response
dopiecewise = 1;

% Cutoff to the penalty for a single outlier
FIXEDLAPSEPDF = 1.5e-6;

% Maximum number of grid refinements in the piecewise linear approximation
MAXITERS = 5;

% Take model parameters
sigmazero = theta(1);
szero = theta(2);
sigmalikezero = theta(3);
slikezero = theta(4);
alpha = theta(5);
beta = theta(6);
sigmaloss = theta(7);
kappa = theta(8);           % Posterior power for stochastic decision making (unused)
lambda = theta(9);
sigmar = theta(13); % Motor/placement noise
wr = theta(14); % Weber's fraction for motor/placement noise

% Number of trials in this condition
nTrials = size(X, 1);

% Increase SSCALE if using MAP estimate
if isinf(kappa) && sigmaloss == 0 && SSCALE < 24; SSCALE = 24; end

% Compute sensory noise std per trial
if szero >= 0
    sigmas = sigmazero*sqrt(1 + (X(:, 2)/szero).^2);
else % Negative szero for cosine noise formula
    % sigmas = sigmazero*(1 + abs(X(:,2)/szero));
    sigmas = sigmazero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(X(:,2)*pi/90))/szero^2);    
end

% Add random jitter to XRANGE to estimate error in the computation of the 
% log likelihood due to the current grid size
% if randomize; xrange = xrange + dx*(rand()-0.5); end

if nargout > 1; doextras = 1; else doextras = 0; end

% Get stimulus grid
srange = getsrange(MAXRNG,SSCALE,sigmalikezero,slikezero,dotaylor);

if ~dopiecewise

    srange_dense = linspace(-MAXRNG,MAXRNG,2000)';
    priorpdf_dense = computePrior(srange_dense,priorinfo);
    priorcdf = cumsum(priorpdf_dense,1);
    priorcdf = priorcdf/priorcdf(end);

    zrange = linspace(0.01,0.99,XGRID);
    [~,idx] = min(abs(bsxfun(@minus,zrange,priorcdf)),[],1);

    % Measurements
    % xrange = alpha*linspace(max(min(X(:,2)-MAXSD*sigmas),-MAXRNG), min(max(X(:,2)+MAXSD*sigmas), MAXRNG), XGRID);
    xrange = srange_dense(idx)';
    xrange = alpha*unique([-MAXRNG:2:xrange(1), xrange, xrange(end):2:MAXRNG]);
    
    % Compute deterministic optimal estimate s*(x)
    if isinf(kappa)
        if dotaylor
            [sstar,extras] = posteriormeantaylor(xrange,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,doextras);
        else
            [sstar,extras] = posteriorestimate(xrange,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,doextras);            
        end
    else
        error('Stochatic decision making not supported.');
    end
    
else
    xrange = []; sstar = [];
    % Starting grid
    xnew = [-MAXRNG:1:-MAXRNG+5,-MAXRNG+10:5:-40,-36:4:-20,-18:2:-10,-9:1:9,10:2:18,20:4:36,40:5:MAXRNG-10,MAXRNG-5:1:MAXRNG];
    TolErr = 0.01;  % Maximum error in piecewise linear approximation
    
    % Loop until estimated error falls below TolErr everywhere
    count = 0;
    while 1
        count = count + 1;
        
        if dotaylor
            [sstarnew,extras] = posteriormeantaylor(xnew,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,0);
        else
            [sstarnew,extras] = posteriorestimate(xnew,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,0);
        end
        [xrange,idx] = sort([xrange, xnew]);
        sstar = [sstar,sstarnew];
        sstar = sstar(idx);
        
        if 1
            hold off;
            scatter(xnew,sstarnew,'ro','FaceColor','r');
            hold on; 
            plot(xrange,sstar,'k','LineWidth',1);
            box off;
            set(gcf,'Color','w');
            set(gca,'TickDir','out');
            xlabel('x (deg)');
            ylabel('s^*(x) (deg)')
            text(-20,20,['Iter #' num2str(count) ', grid size: ' num2str(numel(xrange) + numel(xnew)) '.'],'FontSize',24);
            
            y = computePrior(srange,priorinfo);
            plot(srange,10*y/max(y),'k:','LineWidth',1);
            drawnow;
            
            if any(diff(sstar) < 0)
                theta
                priorinfo
                
                pause;
            end
        end
        
        if count > MAXITERS
            warning('Reached maximum number of grid refinements.');
            break;
        end
        
        % Estimate second derivative
        dx = diff(xrange);        
        h1 = [0,dx];
        h2 = [dx 0];
        f0 = [0 sstar(1:end-1)];        
        f1 = sstar(1:end);
        f2 = [sstar(2:end) 0];
        d2sstar = 2*(h2.*f0 - (h1+h2).*f1 + h1.*f2)./(h1.*h2.*(h1+h2));
        d2sstar([1,end]) = 0;

        % Estimate error bound via Taylor expansion
        err1 = h1.^2.*d2sstar/8;
        err2 = h2.^2.*d2sstar/8;
                
        % Find points where TolErr is violated
        inew = (abs(err1) > TolErr) | (abs(err2) > TolErr);
        xnew = xrange(inew);                
        
        % Break loop if error is less than TolErr everywhere
        if isempty(xnew); break; end
        
        % Add interpolating points
        xnewl = xrange([inew(2:end),false]);
        xnewr = xrange([false,inew(1:end-1)]);        
        % xnew = [xnew + (xnewr-xnew)/3, xnew - (xnew-xnewl)/3];
        xnew = [xnew + 2*(xnewr-xnew)/5, xnew + (xnewr-xnew)/5, xnew - (xnew-xnewl)/5, xnew - 2*(xnew-xnewl)/5];        
    end
end

% numel(xrange)
    

% Piecewise linear approximation of \int p(r|s*(x)) p(x|s) dx
a = xrange(1:end-1);
b = xrange(2:end);
delta = b - a;
alphax = (sstar(2:end) - sstar(1:end-1))./delta;
betax = sstar(1:end-1) - alphax.*b;

%plot(xrange,sstar);
%drawnow

%diff(alphax)
% xrange

% pause

S = X(:, 2)*alpha;              % Remapped stimulus locations
sigmas2 = alpha.^2*sigmas.^2;   % Remapped stimulus noise
R = X(:, 3);                    % Subject's responses

if wr == 0  % Constant motor error
    sigmar2 = sigmar^2;
else        % Weber-like motor error
    sigmar2 = sigmar^2 + wr^2*sstar.^2;
end

% Compute response probability via piecewise linear approximation
sigmasum2 = bsxfun(@plus,sigmar2,bsxfun(@times,alphax.^2,sigmas2));    
zeta = bsxfun_normpdf(bsxfun(@times,alphax,S), bsxfun(@minus,R, betax), sqrt(sigmasum2));
mucdf = bsxfun(@plus, bsxfun(@times,S,sigmar2), bsxfun(@times,bsxfun(@times,alphax,bsxfun(@minus,R,betax)),sigmas2))./sigmasum2;
sigmacdf = sqrt(bsxfun(@rdivide,bsxfun(@times,sigmar2,sigmas2),sigmasum2));
responsepdf = sum(bsxfun(@times,zeta,bsxfun_normcdf(b,mucdf,sigmacdf) - bsxfun_normcdf(a,mucdf,sigmacdf)),2);    

% Add prior-dependent lapse
if lambda > 0
    lapsepdf = computePrior(X(:,3),priorinfo);
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

    %LINSPACE Linearly spaced vector.
    function y = linspace(d1, d2, n)
        y = [d1 + ((0:n-2).*(d2-d1)/(n-1)), d2];
    end

end

%COMPUTEPRIOR Compute unnormalized prior
function y = computePrior(S,priorinfo)

% This is almost zero
NUMZERO = 1e-80;

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
    
%POSTERIORESTIMATE Deterministic estimate based on posterior
function [sstar,extras] = posteriorestimate(xrange,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,doextras)

if isscalar(srange)
    srange = linspace(-MAXRNG, MAXRNG, MAXRNG*srange*2 + 1)';
end
ds = srange(2)-srange(1);

% Compute likelihood
if slikezero >= 0
    sigmasprime = sigmalikezero*sqrt(1 + (srange/slikezero).^2);
else % Negative slikezero for cosine noise formula
    sigmasprime = sigmalikezero.*sqrt(1 + 2*(90/pi)^2*(1 - cos(srange*pi/90))/slikezero^2);
end

% Compute prior        
priorpdf = computePrior(srange,priorinfo);
% priorpdf = priorpdf/(sum(priorpdf, 1)*ds); % Normalize
priorpdf = priorpdf/(simpson1(priorpdf, 1)*ds); % Normalize

% Compute unnormalized posterior for each measurement x in XRANGE
%postpdf = bsxfun(@times, priorpdf, ...
%    bsxfun(@rdivide, exp(-0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2), sigmasprime));


logpostpdf = bsxfun(@plus, log(priorpdf) - log(sigmasprime), ...
    -0.5*bsxfun(@rdivide, bsxfun(@minus, xrange, srange), sigmasprime).^2);
logpostpdf = bsxfun(@minus, logpostpdf, max(logpostpdf,[],1));
postpdf = exp(logpostpdf);

% Deterministic decision making
if sigmaloss == 0                           % MAP
    [~, index] = max(postpdf, [], 1);
    sstar = srange(index)';
elseif sigmaloss == Inf                     % MEAN
    % sstar = sum(bsxfun(@times, srange, postpdf), 1)./sum(postpdf, 1);
    sstar = simpson1(bsxfun(@times, srange, postpdf), 1)./simpson1(postpdf, 1);
elseif sigmaloss == -Inf                    % MEDIAN
    % Linear interpolation of median position from cdf
    cdf = bsxfun(@rdivide,cumtrapz(postpdf, 1),trapz(postpdf, 1));
    [~,pos] = max(cdf >= 0.5, [], 1);
    p0 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + max(pos - 1, 1));
    p1 = cdf(size(cdf,1)*(0:size(cdf,2)-1) + pos);
    coeff1 = (1./(1 + (0.5 - p0)./(p1 - 0.5)))';
    sstar = ((1-coeff1).*srange(pos) + (coeff1).*srange(max(pos-1, 1)))';
end

if doextras
    extras.xrange = xrange;
    extras.srange = srange;
    extras.sstar = sstar;
else
    extras = [];
end        

end


%POSTERIORMEANTAYLOR Deterministic estimate based on posterior
function [sstar,extras] = posteriormeantaylor(xrange,srange,MAXRNG,priorinfo,sigmalikezero,slikezero,sigmaloss,doextras)


% pause

a = srange(1:end-1);
b = srange(2:end);
sbar = 0.5*(a + b);

if priorinfo(3) == 1
    priormu = priorinfo(1);
    priorsigma = priorinfo(2);
    priorw = 1;
else
    priormu(1,1,:) = [priorinfo(1),priorinfo(4)];
    priorsigma(1,1,:) = [priorinfo(2),priorinfo(5)];
    priorw(1,1,:) = [priorinfo(3),1-priorinfo(3)];
end
if priorinfo(6) > 0
    error('Uniform prior not supported with Taylor expansion approximation.');
end

% Compute likelihood
if slikezero >= 0
    sigmas2 = sigmalikezero.^2.*(1 + (sbar/slikezero).^2);
    deltas = 2*(sigmalikezero/slikezero).^2.*sbar;
else % Negative slikezero for cosine noise formula
    sigmas2 = sigmalikezero.^2.*(1 + 2*(90/pi)^2*(1 - cos(sbar*pi/90))/slikezero^2);
    deltas = NaN; 
end

omega = bsxfun(@plus, sigmas2, bsxfun(@times, deltas, bsxfun(@plus, xrange, sbar)));
xi = bsxfun(@times,sigmas2,xrange) + bsxfun(@times, sbar.*deltas, xrange);
tau2 = sigmas2.^3;
zeta = bsxfun(@plus, sigmas2, sbar.*deltas);

m = bsxfun(@rdivide, bsxfun(@plus, bsxfun(@times,priormu,tau2), bsxfun(@times,omega.*xi,priorsigma.^2)), ...
    bsxfun(@plus, tau2, bsxfun(@times,omega.^2,priorsigma.^2)));
V = bsxfun(@rdivide, bsxfun(@times, tau2, priorsigma.^2), ...
    bsxfun(@plus, tau2, bsxfun(@times, omega.^2, priorsigma.^2)));
M = bsxfun_normpdf(xi, bsxfun(@times, omega, priormu), sqrt(bsxfun(@plus, tau2, bsxfun(@times, omega.^2, priorsigma.^2))));
Kzero = 0.5*(erf(bsxfun(@rdivide, bsxfun(@minus, b, m), sqrt(2*V))) - ...
    erf(bsxfun(@rdivide, bsxfun(@minus, a, m), sqrt(2*V))));
Na = bsxfun_normpdf(a,m,sqrt(V));
Nb = bsxfun_normpdf(b,m,sqrt(V));
Izero = sum(bsxfun(@times,priorw, ...
    sum(bsxfun(@times, M, bsxfun(@minus,zeta,bsxfun(@times,deltas,m)).*Kzero ...
    + bsxfun(@times,deltas,bsxfun(@times,V,Nb-Na))),1)),3);
Ione = sum(bsxfun(@times,priorw, ...
    sum(bsxfun(@times, M, bsxfun(@minus, bsxfun(@times,m,zeta), bsxfun(@times,deltas,m.^2+V)).*Kzero + ...
    V.*( bsxfun(@minus, bsxfun(@times, bsxfun(@plus, b, m), deltas), zeta).*Nb - ...
    bsxfun(@minus, bsxfun(@times, bsxfun(@plus, a, m), deltas), zeta).*Na )),1)),3);

sstar = Ione./Izero;

if doextras
    extras.xrange = xrange;
    extras.srange = srange;
    extras.sstar = sstar;
else
    extras = [];
end        

end

%GETSRANGE Get stimulus grid
function srange = getsrange(MAXRNG,srange,sigmalikezero,slikezero,dotaylor)

if ~dotaylor
    if isscalar(srange)
        srange = linspace(-MAXRNG, MAXRNG, MAXRNG*srange*2 + 1)';
    end
    
else

    TolErr = 0.001;   % Tolerance on the approximation error
    w2 = (sigmalikezero/slikezero)^2;
    
    if 0
        srange = 0;
        deltas = [0.1250, 0.1574, 0.1983, 0.2497, 0.3145, 0.3960, 0.4988, 0.6281, 0.7911, 0.9963, 1.2547, 1.5802, 1.9900, 2.5062, 3.1564, 3.9751, 5.0062, 6.3048, 7.9403, 10.0000];
        while 1
            a = srange(end);
            sigmaa = sqrt(sigmalikezero^2 + w2.*a.^2);
            stest = a + deltas;
            sigmas = sqrt(sigmalikezero^2 + w2.*stest.^2);

            err = abs(abs(1./sigmas - (stest - a).*w2.*stest./sigmas.^3 - 1./sigmaa) - TolErr);
            [~,idx] = min(err,[],2);
            idx = min(idx);
            srange(end+1) = stest(idx);
            if srange(end) >= MAXRNG
                srange(end) = MAXRNG;
                srange = [-srange(end:-1:2), srange]';
                break;
            end
        end
        
    else
        
        srange = [0.5,1:1:9,10:2:18,20:4:36,40:5:MAXRNG-10,MAXRNG-5:1:MAXRNG];
        
        count = 0;
        while 1
            count = count + 1;
            a = srange(1:end-1);
            b = srange(2:end);
            stest = 0.5*(a + b);
            sigmaa = sqrt(sigmalikezero^2 + w2.*a.^2);
            sigmas = sqrt(sigmalikezero^2 + w2.*stest.^2);
            err = abs(1./sigmas - (a - stest).*w2.*stest./sigmas.^3 - 1./sigmaa);

            % Find points where TolErr is violated
            inew = (abs(err) > TolErr);
            snew = srange(inew);
                                    
            % Break loop if error is less than TolErr everywhere
            if isempty(snew); break; end

            % Add interpolating points
            snewl = srange([inew(2:end),false]);
            snewr = srange([false,inew(1:end)]);
            snew = [snew + (snewr-snew)/3, snew + 2*(snewr-snew)/3];
            % snew = [snew + 2*(snewr-snew)/5, snew + (snewr-snew)/5, snew - (snew-snewl)/5, snew - 2*(snew-snewl)/5];        
            
            srange = sort([srange,snew]);
            if count > 5; break; end
            
            if 0
                hold off;
                y = 1./sqrt(sigmalikezero^2 + w2.*snew.^2);
                scatter(snew,y,'ro','FaceColor','r');
                hold on;
                y2 = 1./sqrt(sigmalikezero^2 + w2.*srange.^2);
                plot(srange,y2,'k','LineWidth',1);
                box off;
                set(gcf,'Color','w');
                set(gca,'TickDir','out');
                xlabel('s (deg)');
                ylabel('1/\sigma(s) (deg^{-1})')
                text(30,max(y2),['Iter #' num2str(count) ', grid size: ' num2str(numel(srange)) '.'],'FontSize',24);
                drawnow; 
                % pause;
            end
            
            
        end
        srange = [-srange(end:-1:2), srange]';
    end
end
    
end