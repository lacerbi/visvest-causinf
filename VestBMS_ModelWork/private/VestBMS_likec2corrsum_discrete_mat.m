function likec2 = VestBMS_likec2corrsum_discrete(priorpdf2d,like_vis,like_vest)
%VESTBMS_LIKEC2CORRSUM_DISCRETE Compute p(x_vis,x_vest|C=2) for DISCRETE,CORRELATED prior
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-1] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
% SRANGE_VEST: vestibular stimuli. [S-by-1] (double)
%
% ================ OUTPUT VARIABLES ==================
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

likec2(:,:) = sum(bsxfun(@times,bsxfun(@times,priorpdf2d,like_vis),like_vest),1);