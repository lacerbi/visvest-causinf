function postright = VestBMS_PostRight(postpdf)
%VESTBMD_POSTRIGHT Posterior probability of perceiving rightward motion.

n = size(postpdf,1);
pmin = realmin*n/2;

% postpdf = postpdf + realmin;    % Avoid zero probability
idx0deg = (n+1)/2; % Index of 0 deg
postleft = qtrapz(postpdf(1:idx0deg,:,:),1) + pmin;
postright = qtrapz(postpdf(idx0deg:end,:,:),1) + pmin;
postright = postright./(postright + postleft);

end