function logP = VestBMS_logPrior(mp,infostruct)
%PROJECT_LOGPRIOR Model-dependent log prior for a given parameter struct.

%if nargin < 2; infostruct = []; end
%theta = mp.theta;
%model = mp.model;

logP = 0;

% Noninformative uniform prior on all remaining parameters
% uniformpriormask(isnan(uniformpriormask)) = [];    

% All uniform priors are zero (I could always rescale the range)
% mp.logprior = mp.logprior - sum(log(mp.bounds.UB(uniformpriormask) - mp.bounds.LB(uniformpriormask)));       

end