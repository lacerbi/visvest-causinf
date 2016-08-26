function [postright_c1,likec1] = VestBMS_c1postandlikec1sum_discrete(postpdf_c2_uni,like_vis_uni,srange_uni)
%VESTBMS_C1POSTANDLIKEC1QTRAPZ Multiple computations for C=1 (correlated,discrete)
%
% ================ INPUT VARIABLES ====================
% POSTPDF_C2_UNI: p(s) * p(x_vest|s). [S-by-1-by-K] (double)
% LIKE_VIS_UNI: p(x_vis|s). [S-by-K] (double)
% SRANGE_VEST_UNI: s range. [S-by-1] (double)
%
% ================ OUTPUT VARIABLES ==================
% POSTRIGHT_C1: p(right|x_vis,x_vest,C=1). [K-by-K] (double)
% LIKEC1: p(x_vis,x_vest|C=1). [K-by-K] (double)

idx = srange_uni < 0;
temp_left(:,:) = bsxfunandsum(@times,postpdf_c2_uni(idx,:,:),like_vis_uni(idx,:),1,'sum');

idx = srange_uni > 0;
temp_right(:,:) = bsxfunandsum(@times,postpdf_c2_uni(idx,:,:),like_vis_uni(idx,:),1,'sum');

idx = srange_uni == 0;
temp_zero(:,:) = bsxfunandsum(@times,postpdf_c2_uni(idx,:,:),like_vis_uni(idx,:),1,'sum');

likec1 = temp_left + temp_right + temp_zero;
postright_c1 = (temp_right + 0.5*temp_zero) ./ likec1;