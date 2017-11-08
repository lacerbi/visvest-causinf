function prmat = VestBMS_finalqtrapz(xpdf_vis,xpdf_vest,R)
%VESTBMS_FINALQTRAPZ Marginalize response probability over x_vis and x_vest
%
% ================ INPUT VARIABLES ====================
% XPDF_VIS: p(x_vis|s). [S-by-K] (double)
% XPDF_VEST: p(x_vest|s). [S-by-1-by-K] (double)
% R: p(response|x_vis,x_vest). [1-by-K-by-K] (double)
%
% ================ OUTPUT VARIABLES ==================
% RPDF: p(r|s). [S-by-1] (double)

prmat = qtrapz(qtrapz(bsxfun(@times, bsxfun(@times, xpdf_vis, xpdf_vest), R), 2), 3);