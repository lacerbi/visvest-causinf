% VESTBMS_DATALIKE compute minus log likelihood of a dataset.
% 
%   LOGLIKE = VESTBMS_DATALIKE([],MP,INFOSTRUCT) returns log likelihood 
%   of model-parameter struct MP with extra info INFOSTRUCT.
%
%   [LOGLIKE,EXTRAS] = VESTBMS_DATALIKE(...) also returns a series of 
%   additional structures in EXTRAS.
%
%  By Luigi Acerbi <luigi.acerbi@gmail.com>

function [loglike,extras] = VestBMS_datalike(X,mp,infostruct,flags)

if nargin < 4 || isempty(flags); flags = [0 0]; end
if size(flags,2) < 2; flags(2) = 0; end

debug = flags(1);
randomize = flags(2);

extras = [];

% Persistent variables
persistent oldparams;
persistent oldloglikes;
persistent oldtrialloglikes;
persistent starttrials;
persistent ntrials;
if isempty(oldparams)
    % Prepare variables to store parameters between function calls
    for iicnd = 1:4; oldparams{iicnd} = zeros(1, 11); end
    for iicnd = 5:7; oldparams{iicnd} = zeros(1, 26); end
    oldloglikes = zeros(1, 7);
    
    % Count expected number of trial types per condition
    ntrials = zeros(1, 7);
    % for iicnd = 1:4; ntrials(iicnd) = numel(X.unibins{iicnd}(:)); end
    for iicnd = 1:4; ntrials(iicnd) = size(X.unimodal{iicnd},1); end
    for iicnd = 5:7
        for iTask = 1:numel(X.bimbins{iicnd-4})
            % ntrials(iicnd) = ntrials(iicnd) + numel(X.bimbins{iicnd-4}{iTask}(:));
            ntrials(iicnd) = ntrials(iicnd) + size(X.bimodal{iicnd-4}{iTask},1);
        end
    end
    starttrials = [0 cumsum(ntrials)]+1;    
    oldtrialloglikes = NaN(1, sum(ntrials));
end

cnd = mp.cnd;
model = mp.model;

loglikes = zeros(1, max(cnd));
dynamicscale = 1; % Dynamic computation of SSCALE (~ grid points)
% trialloglikes = NaN(1,ntrials(end));

for iicnd = 1:length(cnd)        
    fulltheta = mp.fulltheta{iicnd};   
    
    switch cnd(iicnd)
        case {1, 2, 3, 4} % Unimodal condition
            theta = zeros(1, 8);
            string = {'vis_low', 'vis_med', 'vis_high', 'vest'};
            
            % External noise
            theta(1) = fulltheta.(['sigma_' string{cnd(iicnd)}]);
            theta(2) = fulltheta.(['w_' string{cnd(iicnd)}]);
            if any(cnd(iicnd) == [1 2 3]) && any(model(1) == [3 5 7]); theta(2) = -theta(2);
            elseif cnd(iicnd) == 4 && any(model(2) == [3 5 7]); theta(2) = -theta(2); end
            
            % Internal noise and likelihood
            theta(3) = fulltheta.(['sigmalike_' string{cnd(iicnd)}]);
            theta(4) = fulltheta.(['wlike_' string{cnd(iicnd)}]);
            if any(cnd(iicnd) == [1 2 3]) && any(model(1) == [3 5 7]); theta(4) = -theta(4);
            elseif cnd(iicnd) == 4 && any(model(2) == [3 5 7]); theta(4) = -theta(4); end
            
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
            
            if debug
                [templike extras{cnd(iicnd)}] = VestBMS_UnimodalLeftRightDatalike( ...
                    X.unibins{cnd(iicnd)},model,theta,priorinfo,infostruct.bincenters_uni,infostruct.MAXRNG,mp.XGRID(cnd(iicnd)),SSCALE,0,randomize);
                % templike = log(templike);
                loglikes(cnd(iicnd)) = sum(templike);
                trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = templike;
            else                
                if all(oldparams{cnd(iicnd)} == [theta priorinfo]) 
                    % Do not recompute log likelihood if no changes to local parameter vector
                    loglikes(cnd(iicnd)) = oldloglikes(cnd(iicnd));
                    trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = oldtrialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1);
                else                    
                    DD = X.unimodal{cnd(iicnd)};
                    if ~isempty(DD)
                        templike = VestBMS_UnimodalLeftRightDatalike( ...
                            X.unibins{cnd(iicnd)},model,theta,priorinfo,infostruct.bincenters_uni,infostruct.MAXRNG,mp.XGRID(cnd(iicnd)),SSCALE,0,randomize);
                    else
                        templike = [];
                    end                    
                    % templike = log(templike);
                    loglikes(cnd(iicnd)) = sum(templike);
                    
                    trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = templike;
                end
                
            end
            
        case {5, 6, 7} % Bimodal condition - low, med and high noise
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
            
            maxranges = infostruct.MAXRNG;
            
            % Dynamic assignment of SSCALE
            if dynamicscale
                % [sigmazero_vis,sigmazero_aud,priorinfo(2),priorinfo(5)]
                minsigma = min([theta(1),theta(7),fulltheta.priorsigma])/4;
                SSCALE = min(max(ceil(1/minsigma),1),8);              
            else
                SSCALE = [];
            end
            
            if debug
                if model(9) == 4 || model(9) == 5
                    [templike extras{cnd(iicnd)}] = VestBMS_BimodalLeftRightDatalike_discrete(...
                        X.bimbins{iNoise},model,theta,priorinfo,infostruct.bincenters_bim,mp.XGRID(cnd(iicnd)),0,randomize);                    
                else
                    [templike extras{cnd(iicnd)}] = VestBMS_BimodalLeftRightDatalike(...
                        X.bimbins{iNoise},model,theta,priorinfo,infostruct.bincenters_bim,maxranges,mp.XGRID(cnd(iicnd)),SSCALE,0,randomize);
                end
                % templike = log(templike);
                loglikes(cnd(iicnd)) = sum(templike);
                trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = templike;
            else                
                if all(oldparams{cnd(iicnd)} == [theta priorinfo])
                    % Do not recompute log likelihood if no changes to local
                    % parameter vector
                    loglikes(cnd(iicnd)) = oldloglikes(cnd(iicnd));
                    trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = oldtrialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1);                    
                else
                    if model(9) == 4 || model(9) == 5
                        % Discrete correlated prior
                        templike = VestBMS_BimodalLeftRightDatalike_discrete(...
                            X.bimbins{iNoise},model,theta,priorinfo,infostruct.bincenters_bim,mp.XGRID(cnd(iicnd)),0,randomize);
                    else
                        templike = VestBMS_BimodalLeftRightDatalike(...
                            X.bimbins{iNoise},model,theta,priorinfo,infostruct.bincenters_bim,maxranges,mp.XGRID(cnd(iicnd)),SSCALE,0,randomize);
                    end
                    % templike = log(templike);
                    loglikes(cnd(iicnd)) = sum(templike);                                        
                    trialloglikes(starttrials(cnd(iicnd)):starttrials(cnd(iicnd)+1)-1) = templike;
                end
            end
    end
    
    oldparams{cnd(iicnd)} = [theta priorinfo];
    
end

oldloglikes = loglikes;
oldtrialloglikes = trialloglikes;

if any(isnan(loglikes) | isinf(loglikes) | ~isreal(loglikes))
    loglikes
end

loglike = sum(loglikes);
if debug; temp = extras; extras = []; extras.struct = temp; end
if ~isempty(trialloglikes)
    extras.trialloglikes = trialloglikes(~isnan(trialloglikes));
    extras.ntrials = ntrials;    
end

end