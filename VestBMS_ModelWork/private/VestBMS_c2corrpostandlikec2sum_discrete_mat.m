function [postright_c2,likec2] = VestBMS_c2corrpostandlikec2sum_discrete(priorpdf2d,like_vis,like_vest,srange_vest)
%VESTBMS_C2CORRPOSTANDLIKEC2SUM_DISCRETE Multiple computations for C=2 (correlated,discrete)
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-1] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
% SRANGE_VEST: vestibular stimuli. [S-by-1] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

idx = srange_vest < 0;
temp_left(:,:) = bsxfunandsum(@times,@times,priorpdf2d(idx),like_vis(idx,:),like_vest(idx,:,:),1,'sum');

idx = srange_vest > 0;
temp_right(:,:) = bsxfunandsum(@times,@times,priorpdf2d(idx),like_vis(idx,:),like_vest(idx,:,:),1,'sum');

idx = srange_vest == 0;
temp_zero(:,:) = bsxfunandsum(@times,@times,priorpdf2d(idx),like_vis(idx,:),like_vest(idx,:,:),1,'sum');

likec2 = temp_left + temp_right + temp_zero;
postright_c2 = (temp_right + 0.5*temp_zero) ./ likec2;