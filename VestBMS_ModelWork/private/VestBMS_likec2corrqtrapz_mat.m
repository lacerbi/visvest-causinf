function likec2 = VestBMS_likec2corrqtrapz(priorpdf2d,like_vis,like_vest)
%VESTBMS_LIKEC2CORRQTRAPZ Compute p(x_vis,x_vest|C=2) for CORRELATED prior
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

tmp_vis(:,1,:) = like_vis;
tmp_vest(1,:,1,:) = squeeze(like_vest);
likec2(:,:) = qtrapz(bsxfunandsum(@times,@times,priorpdf2d,tmp_vis,tmp_vest,1,'qtrapz'),2);