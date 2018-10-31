function [theta,LL] = fitDeWinkel2017data
%FITDEWINKEL2017DATA Fit noise distributions to unisensory data

dataMats = load('deWinkel2017data.mat');
data = dataMats.data;
nSubjs = numel(dataMats.data);

LB = [-pi/4,    0.5,0.5,0.5,0.5,    log([1 1 1 1]),         -3,-3,-3,-3,    0];
UB = [pi/4,     2,2,2,2,            log([250 250 250 250]), 1,1,1,1,        0.2];

Nstarts = 101;    % Number of starting points for optimization

for iSubj = 1:nSubjs
    
    % Set starting point from ML solution (different parameterization)
    theta0_orig = data{iSubj}.Params;
    theta0_dewinkel = [theta0_orig([1, 2 5:7, 8 11:13, 14 17:19]), 0.0005];
    theta0_dewinkel(10:13) = theta0_dewinkel(10:13) ./ theta0_dewinkel(6:9);
    theta0_dewinkel(6:9) = log(theta0_dewinkel(6:9));
    
    % Prepare dataset, separate trials by condition (visual and three inertial)
    idx = data{iSubj}.Mat(:,2) == 1;
    X{1} = data{iSubj}.Mat(idx,[4 6]);
    for iCnd = 2:4
        idx = data{iSubj}.Mat(:,2) == 2 & data{iSubj}.Mat(:,3) == (iCnd-1);
        X{iCnd} = data{iSubj}.Mat(idx,[5 6]);
    end
    
    opts.Display = 'iter';
    
    for iModel = 1:2
        theta_best = [];
        fval_best = Inf;
               
        for iStart = 1:Nstarts            
            if iStart == 1
                theta0 = theta0_dewinkel;   % First start from solution from deWinkel et al (2017)
            else
                theta0 = rand(1,numel(LB)).*(UB-LB) + LB;                
            end
            
            fprintf('Fitting subject %d, model %d, restart %d.\n\n', iSubj, iModel, iStart);
            tic
            [theta_iter,fval_iter] = fmincon(@(theta_) fit_nLL(theta_,X,iModel), theta0, [],[],[],[],LB,UB,[],opts);
            toc

            if fval_iter < fval_best
                theta_best = theta_iter;
                fval_best = fval_iter;
            end
        end
        
        theta(iSubj,:,iModel) = theta_best;
        LL(iSubj,iModel) = -fval_best;
        
    end
end

save('noisefits.mat','theta','LL');

end

function nLL = fit_nLL(theta,X,model)
%FIT_NLL Compute negative log likelihood for model fitting

% MODEL is 1 for Gaussian and 2 for von Mises

beta0 = theta(1);                   % Shift bias (common to all conditions)
beta1 = theta(2:5);                 % Scale bias
gamma0 = exp(theta(6:9));           % Base precision/concentration parameter
gamma1 = theta(10:13) .* gamma0;    % Precision/concentration parameter modulation (as per original parameterization)
lambda = theta(14);                 % Lapse rate

% Negative log likelihood
nLL = 0;

% Prepare grid for computing normalization of wrapped Gaussians
if model == 1
    x = linspace(-pi,pi,2e3);
    dx = x(2)-x(1);
end

% Data are already divided by condition
nCnd = numel(X);

for iCnd = 1:nCnd
    theta = X{iCnd}(:,1);   % Stimulus (radians)
    R = X{iCnd}(:,2);       % Response (radians)
    
    % Compute bias and precision/concentration parameter
    mu = beta0 + atan2(beta1(iCnd)*sin(theta),cos(theta));
    kappa = gamma0(iCnd) - gamma1(iCnd)*abs(sin(2*theta));

    switch model
        case 1      % Wrapped Gaussian noise model
            
            sigma = 1./sqrt(kappa);     % Standard deviation from precision

            % Normalization of wrapped normal distributions
            y = bsxfun_normpdf(x,mu,sigma) + bsxfun_normpdf(x,mu-2*pi,sigma) + bsxfun_normpdf(x,mu+2*pi,sigma);
            nf = qtrapz(y,2)*dx;

            % Compute probability of response
            pr = (bsxfun_normpdf(R,mu,sigma) + bsxfun_normpdf(R,mu-2*pi,sigma) + bsxfun_normpdf(R,mu+2*pi,sigma)) ./ nf;
            
        case 2      % von Mises noise model
            
            % Normalize with rescaled Bessel function to avoid numerical errors
            pr = exp(kappa.*cos(R-mu)-kappa) ./ (2*pi*besseli(0,kappa,1));
            
    end
    
    % Add probability of uniform lapsing (takes care of outliers)
    pr_lapse = (1-lambda) * pr + lambda / (2*pi);
    
    nLL = nLL - sum(log(pr_lapse));
end

end