%VESTBMS_C1POSTANDLIKEC1SUM_DISCRETE_TEST
%   Test script for MEX-file VestBMS_c1postandlikec1sum_discrete.
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
srange_uni = 10*rand([S,1]);	%SRANGE_VEST_UNI: s range.

fprintf('============================================\n');
fprintf('Testing VestBMS_c1postandlikec1sum_discrete:\n');
fprintf('============================================\n');

% Call MATLAB and MEX functions
tic; [postright_c1,likec1] = VestBMS_c1postandlikec1sum_discrete_mat(postpdf_c2_uni,like_vis_uni,srange_uni); t = toc;
tic; [postright_c1_mex,likec1_mex] = VestBMS_c1postandlikec1sum_discrete(postpdf_c2_uni,like_vis_uni,srange_uni); t_mex = toc;

% Correctness check
postright_c1_err = sum(abs(postright_c1(:)-postright_c1_mex(:)));
fprintf('Total error (postright_c1): %g\n', postright_c1_err);
if postright_c1_err > TolErr*numel(postright_c1);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in postright_c1.');
end
likec1_err = sum(abs(likec1(:)-likec1_mex(:)));
fprintf('Total error (likec1): %g\n', likec1_err);
if likec1_err > TolErr*numel(likec1);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in likec1.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
