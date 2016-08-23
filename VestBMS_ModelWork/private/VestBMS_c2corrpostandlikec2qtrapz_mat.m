function [postright_c2,likec2] = VestBMS_c2corrpostandlikec2qtrapz(priorpdf2d,like_vis,like_vest)
%VESTBMS_C2CORRPOSTANDLIKEC2QTRAPZ Multiple computations for C=2 (correlated) 
%
% ================ INPUT VARIABLES ====================
% PRIORPDF2D: 2d prior over s_vis, s_vest. [S-by-S] (double)
% LIKE_VIS: visual likelihood. [S-by-K] (double)
% LIKE_VEST: vestibular likelihood. [S-by-1-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C2: p(right|x_vis,x_vest,C=2). [K-by-K] (double)
% LIKEC2: p(x_vis,x_vest|C=2). [K-by-K] (double)

tmp_vis(:,1,:) = like_vis;
tmp_vest(1,:,1,:) = squeeze(like_vest);
postpdf_c2(:,:,:) = bsxfunandsum(@times,@times,priorpdf2d,tmp_vis,tmp_vest,1,'qtrapz');
postright_c2(:,:) = VestBMS_PostRight(postpdf_c2);
likec2(:,:) = qtrapz(postpdf_c2,1);