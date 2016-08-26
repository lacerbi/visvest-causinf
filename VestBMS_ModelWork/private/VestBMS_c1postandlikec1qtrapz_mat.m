function [postright_c1,likec1] = VestBMS_c1postandlikec1qtrapz(postpdf_c2,like_vis)
%VESTBMS_C1POSTANDLIKEC1QTRAPZ Multiple computations for C=1 (uncorrelated)
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS: p(x_vis|s). [S-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C1: p(right|x_vis,x_vest,C=1). [K-by-K] (double)
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

postpdf_c1 = bsxfun(@times,postpdf_c2,like_vis);
postright_c1(:,:) = VestBMS_PostRight(postpdf_c1);
likec1(:,:) = qtrapz(postpdf_c1,1);