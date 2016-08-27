%VESTBMS_LIKEC2CORRSUM_DISCRETE_TEST
%   Test script for MEX-file VestBMS_likec2corrsum_discrete.
%
%   Template MATLAB code generated on 26-Aug-2016 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
K = 401;
S = 88;

% Randomly initialize input variables
% (or write here alternative initializations)
priorpdf2d = 10*rand([S,1]);	%PRIORPDF2D: 2d prior over s_vis, s_vest.
like_vis = 10*rand([S,K]);	%LIKE_VIS: visual likelihood.
like_vest = 10*rand([S,1,K]);	%LIKE_VEST: vestibular likelihood.

fprintf('=======================================\n');
fprintf('Testing VestBMS_likec2corrsum_discrete:\n');
fprintf('=======================================\n');

% Call MATLAB and MEX functions
tic; [likec2] = VestBMS_likec2corrsum_discrete_mat(priorpdf2d,like_vis,like_vest); t = toc;
tic; [likec2_mex] = VestBMS_likec2corrsum_discrete(priorpdf2d,like_vis,like_vest); t_mex = toc;

% Correctness check
likec2_err = sum(abs(likec2(:)-likec2_mex(:)));
fprintf('Total error (likec2): %g\n', likec2_err);
if likec2_err > TolErr*numel(likec2);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in likec2.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
