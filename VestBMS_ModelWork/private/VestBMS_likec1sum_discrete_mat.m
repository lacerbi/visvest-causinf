function likec1 = VestBMS_likec1sum_discrete(postpdf_c2_uni,like_vis_uni)
%VESTBMS_LIKEC1SUM_DISCRETE Multiple computations for C=1 (correlated,discrete)
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2_UNI: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS_UNI: p(x_vis|s). [S-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

likec1(:,:) = sum(bsxfun(@times,postpdf_c2_uni,like_vis_uni),1);