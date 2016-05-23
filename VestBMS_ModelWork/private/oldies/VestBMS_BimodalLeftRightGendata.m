% CUEBMS_BIMODALGENDATA Generate N new datasets from bimodal dataset X
% under model MODEL, parameter set THETA and prior PRIORINFO:
%
% X is a unimodal data matrix. Each row is a trial. For a given row, the 
% columns contain data for:
% X(1) Trial number (unused),
% X(2) Stimulus position (deg),
% X(3) Response (deg) (unused).
%
% MODEL is the model class parameter vector.
%
% For THETA and PRIORINFO, see CUEBMS_UNIMODALDATALIKE.
% 
% Example of external usage:
% d = data{1, 3}; model = [6 3 5 2 1]; theta = [0.03 -0.055, 0.06 0.14, 10 0 0.2];
% ParticleCatch_datalike(d.niceData, d.priormix, model, theta, d.priorsinglegauss)
%
function R = VestBMS_BimodalLeftRightGendata(X,N,model,theta,priorinfo,bincenters,MAXRNG,XGRID,SSCALE)

% Program constants
if nargin < 8; XGRID = 201; end
if nargin < 9; SSCALE = 8; end

nTrialTypes = size(X{2}, 1);
nTrialsPerType = sum(X{2}, 2);

% Generate response probability matrix
[~,extras] = VestBMS_BimodalLeftRightDatalike(X,model,theta,priorinfo,bincenters,MAXRNG,XGRID,SSCALE,[]);

% Take probability of responding LEFT
prmat_left = extras.responsepdf(1:size(extras.responsepdf,1)/2);

R = zeros(nTrialTypes,2,N);
for i = 1:nTrialTypes
    left_r = rand(nTrialsPerType(i),N) < prmat_left(i);
    R(i,1,:) = sum(left_r,1);
    R(i,2,:) = nTrialsPerType(i) - R(i,1,:);
end

end