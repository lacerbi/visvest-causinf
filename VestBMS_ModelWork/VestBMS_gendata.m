function gendata = VestBMS_gendata(N,varargin)
%VESTBMS_GENDATA   Generate a number of fake datasets.
%   GENDATA = VESTBMS_GENDATA(N,MFIT) generates N fake datasets data from 
%   individual model struct MFIT. GENDATA is a cell array of generated data 
%   matrices.
%
%   GENDATA = CUEBMS_GENDATA(N,X,MP,INFOSTRUCT) generates fake datasets 
%   from data matrix X, model parameter structure MP and INFOSTRUCT.
%
%  By Luigi Acerbi <luigi.acerbi@gmail.com>

if nargin < 4
    mfit = ModelWork_loadFields('VestBMS',varargin{1});
    X = mfit.X;
    mp = mfit.mp;
    infostruct = mfit.infostruct;
else
    X = varargin{1};
    mp = varargin{2};
    infostruct = varargin{3};
end

gendata = [];
if ~isempty(mfit.sampling) && ~isempty(mfit.sampling.samples)
    Nsamples = size(mfit.sampling.samples,1);
    idx = round(linspace(1,Nsamples,N));
    for i = 1:N
        mp = VestBMS_setupModel(mp, mfit.sampling.samples(idx(i),:), mfit.model, infostruct);
        XX = GenDatasets(1,mp,X,infostruct);
        gendata{i} = squeeze(XX(:, :, 1));
    end
else
    mp = VestBMS_setupModel(mp, mfit.maptheta, mfit.model, infostruct);    
    XX = GenDatasets(N,mp,X,infostruct);    
    for i = 1:N; gendata{i} = squeeze(XX(:, :, i)); end    
end


end

function XX = GenDatasets(N,mp,X,infostruct)
%GENDATASETS Generate fake datasets.

cnd = mp.cnd;
model = mp.model;

dynamicscale = 1; % Dynamic computation of SSCALE (~ grid points)

XX = [];
clear functions;

for iicnd = 1:length(cnd)
        
    fulltheta = mp.fulltheta{iicnd}; % This is the MAP (mean for sampling)
    
    switch cnd(iicnd)
        case {1, 2, 3, 4} % Unimodal condition
            iindex = cnd(iicnd);
            
            theta = zeros(1, 8);
            string = {'vis_low', 'vis_med', 'vis_high', 'vest'};
            
            % External noise
            theta(1) = fulltheta.(['sigma_' string{iindex}]);
            theta(2) = fulltheta.(['w_' string{iindex}]);
            if any(iindex == [1 2 3]) && any(model(1) == [3 5 7]); theta(2) = -theta(2);
            elseif iindex == 4 && any(model(2) == [3 5 7]); theta(2) = -theta(2); end
            
            % Internal noise and likelihood
            theta(3) = fulltheta.(['sigmalike_' string{iindex}]);
            theta(4) = fulltheta.(['wlike_' string{iindex}]);
            if any(iindex == [1 2 3]) && any(model(1) == [3 5 7]); theta(4) = -theta(4);
            elseif iindex == 4 && any(model(2) == [3 5 7]); theta(4) = -theta(4); end
            
            % Decision making, lapse and errors
            theta(8) = 1/fulltheta.tau_softmax; % Inverse temperature
            
            priorinfo = [fulltheta.priormu,fulltheta.priorsigma];
            
            % Dynamic assignment of SSCALE
            if dynamicscale
                minsigma = min([theta(1),fulltheta.priorsigma])/4;
                SSCALE = min(max(ceil(1/minsigma),1),8);
            else
                SSCALE = [];
            end
            
            % Fill data matrix
            datafun = @(X) VestBMS_UnimodalLeftRightDatalike(X,model,theta,priorinfo,infostruct.bincenters_uni,infostruct.MAXRNG,[],SSCALE);
            Rmat = LeftRightGendata(datafun,X.unibins{iindex},N);
            
            % Convert from response matrix to trial matrix
            R = [];
            for iMat = 1:N
                mat = Rmat(:,:,iMat);
                temp = [];
                for t = 1:size(mat,1)
                    nTot = sum(mat(t,:),2);
                    resp = [-ones(mat(t,1),1); ones(mat(t,2),1)];
                    D = [zeros(nTot,1), infostruct.bincenters_uni(t)*ones(nTot,1), resp];
                    temp = [temp; D];           
                end
                R(:,:,iMat) = temp;
            end
            
            D = zeros(size(R, 1), 10, size(R, 3));
            D(:, 1, :) = R(:, 1, :) + 600; % Trial number
            D(:, 3, :) = mod(iindex, 4); % Visual noise level
            D(:, 7, :) = 1; % Number of stimuli (1 Unimodal, 2 Bimodal)
            if iindex == 4 % Unimodal vestibular
                D(:, 4, :) = 2; % Trial response type (2 Vestibular)
                D(:, 5, :) = R(:, 2, :); % Vestibular stimulus position
                D(:, 8, :) = R(:, 3, :); % Vestibular response
            else % Unimodal video
                D(:, 4, :) = 1; % Trial response type (1 Visual)
                D(:, 6, :) = R(:, 2, :); % Visual stimulus position
                D(:, 9, :) = R(:, 3, :); % Visual response
            end

        case {5, 6, 7} % Bimodal condition
            iNoise = cnd(iicnd) - 4;
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
            end
            if isfield(fulltheta,'tau_causinf')
                theta(16) = fulltheta.tau_causinf;
            end
                                    
            priorinfo = [fulltheta.priormu fulltheta.priorsigma fulltheta.pcommon fulltheta.kcommon];
                        
            % Dynamic assignment of SSCALE
            if dynamicscale
                % [sigmazero_vis,sigmazero_aud,priorinfo(2),priorinfo(5)]
                minsigma = min([theta(1),theta(7),fulltheta.priorsigma])/4;
                SSCALE = min(max(ceil(1/minsigma),1),8);                
            else
                SSCALE = [];
            end
            
            % Fill data matrix
            datafun = @(X) VestBMS_BimodalLeftRightDatalike(X,model,theta,priorinfo,infostruct.bincenters_bim,infostruct.MAXRNG,[],SSCALE);
            Rmat = LeftRightGendata(datafun,X.bimbins{iNoise},N);
            
            % Rmat = VestBMS_BimodalLeftRightGendata(X.bimbins{iNoise},N,model,theta,priorinfo,infostruct.bincenters_bim,infostruct.MAXRNG,SSCALE);
                        
            % Convert from response matrix to trial matrix
            R = [];
            for iMat = 1:N
                mat = Rmat(:,:,iMat);
                temp = [];
                for t = 1:size(mat,1)
                    nTot = sum(mat(t,:),2);
                    resp = [-ones(mat(t,1),1); ones(mat(t,2),1)];
                    D = [zeros(nTot,1), 2*ones(nTot,1), infostruct.bincenters_bim{2}(t)*ones(nTot,1), infostruct.bincenters_bim{1}(t)*ones(nTot,1), resp];
                    temp = [temp; D];           
                end
                R(:,:,iMat) = temp;
            end
            
            D = zeros(size(R, 1), 10, size(R, 3));
            D(:, 1, :) = R(:, 1, :) + 600; % Trial number
            D(:, 3, :) = iNoise; % Visual noise level
            D(:, 4, :) = R(:, 2, :); % Trial response type
            D(:, 5, :) = R(:, 3, :); % Vestibular stimulus position
            D(:, 6, :) = R(:, 4, :); % Visual stimulus position
            D(:, 7, :) = 2; % Number of stimuli (1 Unimodal, 2 Bimodal)
            fvest = D(:, 4, 1) == 2;
            D(fvest, 8, :) = R(fvest, 5, :); % Vestibular response
            fvis = D(:, 4, 1) == 1;
            D(fvis, 9, :) = R(fvis, 5, :); % Visual response
            fcat = D(:, 4, 1) == 3;
            temp = R(fcat, 5, :); % Categorical response
            temp(temp == 0) = 2; % Convert zeros to 2
            D(fcat, 10, :) = temp;
    end
    
    XX = [XX; D];
end

end

%--------------------------------------------------------------------------
function Rmat = LeftRightGendata(datafun,X,N)
%LEFTRIGHTGENDATA Generate left/right unisensory and bisensory datasets.

if iscell(X)
    nTrialTypes = size(X{2}, 1);
    nTrialsPerType = sum(X{2}, 2);
else
    nTrialTypes = size(X, 1);
    nTrialsPerType = sum(X, 2);    
end

% Generate response probability matrix
[~,extras] = datafun(X);

% Take probability of responding LEFT
prmat_left = extras.responsepdf(1:size(extras.responsepdf,1)/2);

Rmat = zeros(nTrialTypes,2,N);
for i = 1:nTrialTypes
    left_r = rand(nTrialsPerType(i),N) < prmat_left(i);
    Rmat(i,1,:) = sum(left_r,1);
    Rmat(i,2,:) = nTrialsPerType(i) - Rmat(i,1,:);
end

end
