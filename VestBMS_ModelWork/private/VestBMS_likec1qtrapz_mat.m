function likec1 = VestBMS_likec1qtrapz(postpdf_c2,like_vis)
%VESTBMS_LIKEC1QTRAPZ Compute p(x_vis,x_vest|C=1) for uncorrelated prior
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS: p(x_vis|s). [S-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

likec1(:,:) = qtrapz(bsxfun(@times,postpdf_c2,like_vis),1);