function [infodata,ll_est] = VestBMS_loglikeCheck(mfit,varargin)
%VESTBMS_LOGLIKECHECK Check validity of computation and data generation.
%   VESTBMS_LOGLIKECHECK(MFIT) tests the computation of the log likelihood
%   and fake data generation for provided model structure MFIT.
%
%   VESTBMS_LOGLIKECHECK(MFIT,N) uses N function evaluations (default N=1e3).
%
%
%   See also MODELWORK_LOGLIKECHECK.

if nargin < 1
    help VestBMS_loglikeCheck;
    return;
end

if nargin < 2 || ...
        ( nargin == 2 && isnumeric(varargin{1}) && isscalar(varargin{1}) && isempty(varargin{1}))
    if nargin < 2 || isempty(varargin{1}); N = []; else N = varargin{1}; end
    ModelWork_loglikeCheck('VestBMS',mfit,N);
    return;
else
    mfit_true = varargin{1};
    infodata = varargin{2};
end

lambda = 0.0001;
X = mfit.X;

if isempty(infodata)
    infodata.Rmats = [];
    firstcall = 1;
else
    firstcall = 0;
end
idx = size(infodata.Rmats,4) + 1;

if ~isempty(X.bimbins{1}{2}); task = 2; else task = 3; end
for iNoise = 1:3
    infodata.Rmats(:,:,iNoise,idx) = X.bimbins{iNoise}{task};
end

if firstcall
    % True response matrix from data
    for iNoise = 1:3
        infodata.trueRmat(:,:,iNoise) = mfit_true.X.bimbins{iNoise}{task};
    end
end

% Approximate likelihood via effective samples
Rmean = mean(infodata.Rmats(:,:,:,1:idx),4);
Rmean = bsxfun(@rdivide, Rmean, Rmean(:,1,:)+Rmean(:,2,:));
Rmean = lambda*0.5 + (1-lambda)*Rmean;

ll_est = nansum(log(Rmean(:)).*infodata.trueRmat(:));

end