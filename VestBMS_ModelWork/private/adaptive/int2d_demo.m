% CUEBMS_BIMODALDATALIKE Calculate (log) likelihood of bimodal dataset X
function int2d_demo(mfit,iNoise)

if nargin < 2 || isempty(iNoise); iNoise = 1; end

TolErr = 0.001; % Relative error tolerance

dynamicscale = 1;   % Dynamic scale for SSRANGE
XGRID = 401;

%% DATALIKE parameter initialization
fulltheta = mfit.mp.fulltheta{iNoise};
X = mfit.X.bimbins{iNoise};
model = mfit.model;
infostruct = mfit.infostruct;
bincenters = infostruct.bincenters_bim;

theta = zeros(1, 16);
string = {'vis_low', 'vis_med', 'vis_high'};

% External visual noise
theta(1) = fulltheta.(['sigma_' string{iNoise}]);
theta(2) = fulltheta.(['w_' string{iNoise}]);
if any(model(1) == [3 5 7]); theta(2) = -theta(2); end

% Internal visual noise and likelihood
theta(3) = fulltheta.(['sigmalike_' string{iNoise}]);
theta(4) = fulltheta.(['wlike_' string{iNoise}]);
if any(model(1) == [3 5 7]); theta(4) = -theta(4); end

% External vestibular noise
theta(7) = fulltheta.sigma_vest;
theta(8) = fulltheta.w_vest;
if any(model(2) == [3 5 7]); theta(8) = -theta(8); end

% Internal vestibular noise and likelihood
theta(9) = fulltheta.sigmalike_vest;
theta(10) = fulltheta.wlike_vest;
if any(model(2) == [3 5 7]); theta(10) = -theta(10); end

% Decision making, lapse and errors
theta(14) = 1/fulltheta.tau_softmax;

if isfield(fulltheta,'invgamma_causinf')
    theta(15) = 1./fulltheta.invgamma_causinf;
    if ~isfield(fulltheta,'invgamma_causinf_unity')
        fulltheta.invgamma_causinf_unity = fulltheta.invgamma_causinf;
    end
    theta(17) = 1./fulltheta.invgamma_causinf_unity;
end
if isfield(fulltheta,'tau_causinf')
    theta(16) = fulltheta.tau_causinf;
end
theta(18) = fulltheta.lambda;

priorinfo = [fulltheta.priormu fulltheta.priorsigma fulltheta.priorsigmadelta ...
    fulltheta.pcommon fulltheta.kcommon fulltheta.pcommon_unity fulltheta.kcommon_unity];

maxranges = infostruct.MAXRNG;

% Dynamic assignment of SSCALE
if dynamicscale
    % [sigmazero_vis,sigmazero_aud,priorinfo(2),priorinfo(5)]
    minsigma = min([theta(1),theta(7),fulltheta.priorsigma])/4;
    SSCALE = min(max(ceil(1/minsigma),1),8)      
else
    SSCALE = 1;
end

SSCALE = 8;

%% BIMODALLEFTRIGHTDATALIKE computations

% When integrating a Gaussian, go up to this SDs away
MAXSD = 5;

% Take model parameters
sigmazero_vis = theta(1);
w_vis = theta(2);
sigmalikezero_vis = theta(3);
wlike_vis = theta(4);

sigmazero_vest = theta(7);
w_vest = theta(8);
sigmalikezero_vest = theta(9);
wlike_vest = theta(10);

priorsigmadelta = priorinfo(3);

MAXRNG = maxranges(1);
MAXRNG_XMEAS = 180;

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
xrange_vis(1, :, 1) = linspace(max(min(bincenters_vis-MAXSD*sigmas_vis),-MAXRNG_XMEAS), min(max(bincenters_vis+MAXSD*sigmas_vis), MAXRNG_XMEAS), XGRID);
xrange_vest (1, 1, :) = linspace(max(min(bincenters_vest-MAXSD*sigmas_vest),-MAXRNG_XMEAS), min(max(bincenters_vest+MAXSD*sigmas_vest), MAXRNG_XMEAS), XGRID);

% Wrap large noisy measurement around circle
wraparound = MAXRNG_XMEAS >= 180 && ...
        ( min(bincenters_vis-MAXSD*sigmas_vis) <= -180 || max(bincenters_vis+MAXSD*sigmas_vis) >= 180 || ...
        min(bincenters_vest-MAXSD*sigmas_vest) <= -180 || max(bincenters_vest+MAXSD*sigmas_vest) >= 180);

srange = linspace(-MAXRNG, MAXRNG, 1e3)';
ds = diff(srange(1:2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood; variable range depends on type

likerange_vis = srange; 
likerange_vest = srange; 

% Compute sensory likelihood std for vision
if wlike_vis >= 0; likemodel_vis = 'A'; else likemodel_vis = 'C'; end
sigmasprime_vis = VestBMS_sensoryNoise(likemodel_vis,likerange_vis,sigmalikezero_vis,wlike_vis);

% Compute sensory likelihood std for vestibular
if wlike_vest >= 0; likemodel_vest = 'A'; else likemodel_vest = 'C'; end
sigmasprime_vest = VestBMS_sensoryNoise(likemodel_vest,likerange_vest,sigmalikezero_vest,wlike_vest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute likelihood for non-Gaussian likelihoods

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

% Consider only subset of XGRID
idx = 1:50:XGRID;
like_vis = like_vis(:,idx);
like_vest = like_vest(:,1,idx);

% Compute CORRELATED prior, p(s_vis, s_vest), for eccentric noise
srange_vis = srange(:);
srange_vest = srange(:)';
priorpdf2d = bsxfun_normpdf(0.5*srange_vest, -0.5*srange_vis + priorinfo(1), priorinfo(2)) .* ...
    bsxfun_normpdf(srange_vest, srange_vis, priorsigmadelta);
priorpdf2d = priorpdf2d/(qtrapz(qtrapz(priorpdf2d,1))*ds*ds); % Normalize prior

like_vis_res(:,1,:,1) = like_vis;
like_vest_res(1,:,1,:) = like_vest;
int2d = bsxfun(@times,bsxfun(@times,priorpdf2d,like_vis_res),like_vest_res);

tic
% [postright_c2(1,:,:),likec2(1,:,:)] = VestBMS_c2corrpostandlikec2sum(priorpdf2d,like_vis_res,like_vest_res);
[postright_c2(1,:,:),likec2(1,:,:)] = VestBMS_c2corrpostandlikec2qtrapzhere(priorpdf2d,like_vis_res,like_vest_res);
% likec2 = likec2*ds*ds + realmin;
toc

tic
like_vis([1 end],:) = 0.5*like_vis([1 end],:);
like_vest([1 end],:) = 0.5*like_vest([1 end],:);
%[postright_c2bis,likec2_bis,err,fevals] = VestBMS_c2corrpostandlikec2adapt(priorpdf2d,like_vis,like_vest,TolErr);
toc

tic
[postright_c2bis_mex,likec2_bis_mex,err_mex,fevals_mex] = VestBMS_c2corrpostandlikec2adapt_mex(priorpdf2d,like_vis,like_vest);
toc

tic
[postright_c2bis_mex2,likec2_bis_mex2,err_mex2,fevals_mex2] = VestBMS_c2corrpostandlikec2adapt_mex2(priorpdf2d,like_vis,like_vest,TolErr);
toc

squeeze(likec2) ./ likec2_bis

end


function [postright_c2,likec2] = VestBMS_c2corrpostandlikec2sum(priorpdf2d,like_vis,like_vest)
%VESTBMS_C2CORRPOSTANDLIKEC2QTRAPZ Multiple computations for C=2 (correlated) 
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

postpdf_c2(:,:,:) = bsxfunandsum(@times,@times,priorpdf2d,like_vis,like_vest,1,'sum');
postright_c2(:,:) = VestBMS_PostRight(postpdf_c2);
likec2(:,:) = sum(postpdf_c2,1);

end


function [postright_c2,likec2] = VestBMS_c2corrpostandlikec2qtrapzhere(priorpdf2d,like_vis,like_vest)
%VESTBMS_C2CORRPOSTANDLIKEC2QTRAPZ Multiple computations for C=2 (correlated) 
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

postpdf_c2(:,:,:) = bsxfunandsum(@times,@times,priorpdf2d,like_vis,like_vest,1,'qtrapz');

n = size(postpdf_c2,1);
pmin = realmin*n/2;
idx0deg = n/2; % Index of 0 deg

% Use MEX files
postleft = sum(postpdf_c2(1:idx0deg,:,:),1) - 0.5*postpdf_c2(1,:,:) + pmin;
posttemp = sum(postpdf_c2(idx0deg+1:end,:,:),1) - 0.5*postpdf_c2(end,:,:) + pmin;
postright_c2(:,:) = posttemp./(posttemp + postleft);

likec2(:,:) = qtrapz(postpdf_c2,1);

end