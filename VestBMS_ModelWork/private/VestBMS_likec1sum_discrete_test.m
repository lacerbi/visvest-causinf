%VESTBMS_LIKEC1SUM_DISCRETE_TEST
%   Test script for MEX-file VestBMS_likec1sum_discrete.
%
%   Template MATLAB code generated on 26-Aug-2016 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
K = 100;
S = 150;

% Randomly initialize input variables
% (or write here alternative initializations)
postpdf_c2_uni = 10*rand([S,1,K]);	%POSTPDF_C2_UNI: p(s) * p(x_vest|s).
like_vis_uni = 10*rand([S,K]);	%LIKE_VIS_UNI: p(x_vis|s).

fprintf('===================================\n');
fprintf('Testing VestBMS_likec1sum_discrete:\n');
fprintf('===================================\n');

% Call MATLAB and MEX functions
tic; [likec1] = VestBMS_likec1sum_discrete_mat(postpdf_c2_uni,like_vis_uni); t = toc;
tic; [likec1_mex] = VestBMS_likec1sum_discrete(postpdf_c2_uni,like_vis_uni); t_mex = toc;

% Correctness check
likec1_err = sum(abs(likec1(:)-likec1_mex(:)));
fprintf('Total error (likec1): %g\n', likec1_err);
if likec1_err > TolErr*numel(likec1);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in likec1.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
