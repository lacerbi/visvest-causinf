function postright = VestBMS_PostRight(postpdf)
%VESTBMS_POSTRIGHT Posterior probability of perceiving rightward motion.

n = size(postpdf,1);
pmin = realmin*n/2;

idx0deg = (n+1)/2; % Index of 0 deg

% Use MEX files
postleft = qtrapzc(postpdf,1,[1;idx0deg]) + pmin;
posttemp = qtrapzc(postpdf,1,[idx0deg;size(postpdf,1)]) + pmin;
postright(1,:,:) = posttemp./(posttemp + postleft);
    
    %t1 = postright;
    %clear postright;

    %postleft = qtrapz(postpdf(1:idx0deg,:,:),1) + pmin;
    %postright = qtrapz(postpdf(idx0deg:end,:,:),1) + pmin;
    %postright = postright./(postright + postleft);
    
    %sum(abs(t1(:) - postright(:)))

end