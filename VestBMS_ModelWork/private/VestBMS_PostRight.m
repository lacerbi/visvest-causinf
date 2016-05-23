function postright = VestBMS_PostRight(postpdf)
%VESTBMD_POSTRIGHT Posterior probability of perceiving rightward motion.

postpdf = postpdf + realmin;    % Avoid zero probability
idx0deg = (size(postpdf,1)+1)/2; % Index of 0 deg
postleft = qtrapz(postpdf(1:idx0deg,:,:),1);
postright = qtrapz(postpdf(idx0deg:end,:,:),1);
postright = postright./(postright + postleft);

end