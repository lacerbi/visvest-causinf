function XX = VestBMS_gendata(N,mp,X,infostruct,regenerate)
%VESTBMS_GENDATA   Generate a number of fake datasets.
%   XX = VESTBMS_GENDATA(N,MP,X,INFOSTRUCT) generates N fake datasets data 
%   with model parameter structure MP, original dataset X and extra info
%   struct INFOSTRUCT.
%
%   XX = VESTBMS_GENDATA(N,MP,X,INFOSTRUCT,1) regenerates the REAL dataset 
%   (trial order might be swapped).
%
%   See also MODELWORK_GENDATA.
 
cnd = mp.cnd;
model = mp.model;

dynamicscale = 1; % Dynamic computation of SSCALE (~ grid points)

XX = [];

for iicnd = 1:length(cnd)
    clear functions;
        
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
            theta(9) = fulltheta.lambda;
            
            priorinfo = [fulltheta.priormu,fulltheta.priorsigma];
            
            % Dynamic assignment of SSCALE
            if dynamicscale
                minsigma = min([theta(1),fulltheta.priorsigma])/4;
                SSCALE = min(max(ceil(1/minsigma),1),8);
            else
                SSCALE = [];
            end
            
            % Fill data matrix
            datafun = @(X) VestBMS_UnimodalLeftRightDatalike(X,model,theta,priorinfo,infostruct.bincenters_uni,infostruct.MAXRNG,[],SSCALE,[1,N]);
            if regenerate
                Rmat = X.unibins{iindex};
            else
                Rmat = GendataFun(datafun,X.unibins{iindex},N);
            end
            
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
            theta = zeros(1, 18);
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
            
            % Retrocompatibility
            if ~isfield(fulltheta,'pcommon_unity'); fulltheta.pcommon_unity = fulltheta.pcommon; end
            if ~isfield(fulltheta,'kcommon_unity'); fulltheta.kcommon_unity = fulltheta.kcommon; end
            if ~isfield(fulltheta,'priorsigmadelta'); fulltheta.priorsigmadelta = 0; end
            if ~isfield(fulltheta,'priorsigmadelta_unity'); fulltheta.priorsigmadelta_unity = fulltheta.priorsigmadelta; end

            priorinfo = [fulltheta.priormu fulltheta.priorsigma fulltheta.priorsigmadelta fulltheta.priorsigmadelta_unity...
                fulltheta.pcommon fulltheta.kcommon fulltheta.pcommon_unity fulltheta.kcommon_unity];
            
            % Random unity judgments
            if model(16) == 4
                string = {'low', 'med', 'high'};                
                theta(17) = fulltheta.(['random_unity_' string{iNoise}]);             
            end
            
            % Dynamic assignment of SSCALE
            if dynamicscale
                % [sigmazero_vis,sigmazero_aud,priorinfo(2),priorinfo(5)]
                minsigma = min([theta(1),theta(7),fulltheta.priorsigma])/4;
                SSCALE = min(max(ceil(1/minsigma),1),8);                
            else
                SSCALE = [];
            end
            
            % Fill data matrix
            if model(9) == 4 || model(9) == 5
                datafun = @(X) VestBMS_BimodalLeftRightDatalike_discrete(X,model,theta,priorinfo,infostruct.bincenters_bim,[],[1,N]);
            else
                datafun = @(X) VestBMS_BimodalLeftRightDatalike(X,model,theta,priorinfo,infostruct.bincenters_bim,infostruct.MAXRNG,[],SSCALE,[1,N]);
            end
            if regenerate
                Rmat = X.bimbins{iNoise}{2};
                Rmat_unity = X.bimbins{iNoise}{3};                
            else
                [Rmat,Rmat_unity] = GendataFun(datafun,X.bimbins{iNoise},N);
            end
            
            % Rmat = VestBMS_BimodalLeftRightGendata(X.bimbins{iNoise},N,model,theta,priorinfo,infostruct.bincenters_bim,infostruct.MAXRNG,SSCALE);
            
            % Left/right localization trials
            if ~isempty(Rmat)
                R = [];
                % Convert from response matrix to trial matrix
                for iMat = 1:N
                    mat = Rmat(:,:,iMat)';
                    if 0
                        temp = [];
                        for t = 1:size(mat,1)
                            nTot = sum(mat(t,:),2);
                            resp = [-ones(mat(t,1),1); ones(mat(t,2),1)];
                            D = [zeros(nTot,1), 2*ones(nTot,1), infostruct.bincenters_bim{2}(t)*ones(nTot,1), infostruct.bincenters_bim{1}(t)*ones(nTot,1), resp];
                            temp = [temp; D];           
                        end
                    else
                        nTot = sum(mat,1);
                        nTotal = sum(nTot);
                        Rtemp = zeros(nTotal,5);
                        Rtemp(:,2) = 2;     % Task
                        Rtemp(:,3) = repelem(infostruct.bincenters_bim{2},nTot);
                        Rtemp(:,4) = repelem(infostruct.bincenters_bim{1},nTot);
                        Rlist = [-ones(1,size(mat,2)); ones(1,size(mat,2))];
                        Rtemp(:,5) = repelem(Rlist(:),mat(:));
                        R(:,:,iMat) = Rtemp;
                    end
                    R(:,:,iMat) = Rtemp;
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
            else
                D = [];
            end
            
            % Unity judgements trials
            if ~isempty(Rmat_unity)
                R = [];
                % Convert from response matrix to trial matrix
                for iMat = 1:N
                    mat = Rmat_unity(:,:,iMat)';
                    if 0
                        temp = [];
                        for t = 1:size(mat,1)
                            nTot = sum(mat(t,:),2);
                            resp = [ones(mat(t,1),1); zeros(mat(t,2),1)];
                            Du = [zeros(nTot,1), 3*ones(nTot,1), infostruct.bincenters_bim{2}(t)*ones(nTot,1), infostruct.bincenters_bim{1}(t)*ones(nTot,1), resp];
                            temp = [temp; Du];           
                        end
                        R(:,:,iMat) = temp;
                    else
                        nTot = sum(mat,1);
                        nTotal = sum(nTot);
                        Rtemp = zeros(nTotal,5);
                        Rtemp(:,2) = 3;     % Task
                        Rtemp(:,3) = repelem(infostruct.bincenters_bim{2},nTot);
                        Rtemp(:,4) = repelem(infostruct.bincenters_bim{1},nTot);
                        Rlist = [ones(1,size(mat,2)); zeros(1,size(mat,2))];
                        Rtemp(:,5) = repelem(Rlist(:),mat(:));                        
                        R(:,:,iMat) = Rtemp;
                    end
                end

                Du = zeros(size(R, 1), 10, size(R, 3));
                Du(:, 1, :) = R(:, 1, :) + 600; % Trial number
                Du(:, 3, :) = iNoise; % Visual noise level
                Du(:, 4, :) = R(:, 2, :); % Trial response type
                Du(:, 5, :) = R(:, 3, :); % Vestibular stimulus position
                Du(:, 6, :) = R(:, 4, :); % Visual stimulus position
                Du(:, 7, :) = 2; % Number of stimuli (1 Unimodal, 2 Bimodal)
                fvest = Du(:, 4, 1) == 2;
                Du(fvest, 8, :) = R(fvest, 5, :); % Vestibular response
                fvis = Du(:, 4, 1) == 1;
                Du(fvis, 9, :) = R(fvis, 5, :); % Visual response
                fcat = Du(:, 4, 1) == 3;
                temp = R(fcat, 5, :); % Categorical response
                temp(temp == 0) = 2; % Convert zeros to 2
                Du(fcat, 10, :) = temp;
            else
                Du = [];
            end
            
            D = [D; Du];
            
    end
    
    XX = [XX; D];
end

end

%--------------------------------------------------------------------------
function [Rmat,Rmat_unity] = GendataFun(datafun,X,N)
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

% Compute left/right localization
if ~isempty(extras.responsepdf)

    if ~isfield(extras,'Rmat') || isempty(extras.Rmat)    
        % Take probability of responding LEFT
        prmat_left = extras.responsepdf(:,1);

        Rmat = zeros(nTrialTypes,2,N);
        for i = 1:nTrialTypes
            left_r = rand(nTrialsPerType(i),N) < prmat_left(i);
            Rmat(i,1,:) = sum(left_r,1);
            Rmat(i,2,:) = nTrialsPerType(i) - Rmat(i,1,:);
        end
    else
        Rmat = extras.Rmat;
    end
        
else
    Rmat = [];
end

% Compute unity judgements
if nargout > 1 && ~isempty(extras.responsepdf_unity)
    if iscell(X)
        nTrialTypes = size(X{3}, 1);
        nTrialsPerType = sum(X{3}, 2);
    else
        nTrialTypes = size(X, 1);
        nTrialsPerType = sum(X, 2);    
    end
    
    prmat_unity = extras.responsepdf_unity(:,1);

    if ~isfield(extras,'Rmat_unity') || isempty(extras.Rmat_unity)    
        Rmat_unity = zeros(nTrialTypes,2,N);
        for i = 1:nTrialTypes
            unity_r = rand(nTrialsPerType(i),N) < prmat_unity(i);
            Rmat_unity(i,1,:) = sum(unity_r,1);
            Rmat_unity(i,2,:) = nTrialsPerType(i) - Rmat_unity(i,1,:);
        end
    else
        Rmat_unity = extras.Rmat_unity;
    end
else
    Rmat_unity = [];
end


end
