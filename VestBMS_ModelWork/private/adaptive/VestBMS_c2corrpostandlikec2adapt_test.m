%VESTBMS_C2CORRPOSTANDLIKEC2ADAPT_TEST
%   Test script for MEX-file VestBMS_c2corrpostandlikec2adapt.
%
%   Template MATLAB code generated on 25-Aug-2016 with MEXXER v0.2 
%   (https://github.com/lacerbi/mexxer).

TolErr = sqrt(eps);	% Maximum error tolerance per array element

% Define array sizes (choose reasonable values)
K = 100;
S = 150;

% Randomly initialize input variables
% (or write here alternative initializations)
priorpdf2d = 10*rand([S,S]);	%PRIORPDF2D: 2d prior over s_vis, s_vest.
like_vis = 10*rand([S,K]);	%LIKE_VIS: visual likelihood.
like_vest = 10*rand([S,1,K]);	%LIKE_VEST: vestibular likelihood.

fprintf('=========================================\n');
fprintf('Testing VestBMS_c2corrpostandlikec2adapt:\n');
fprintf('=========================================\n');

% Call MATLAB and MEX functions
tic; [postright_c2,likec2,err,fevals] = VestBMS_c2corrpostandlikec2adapt(priorpdf2d,like_vis,like_vest); t = toc;
tic; [postright_c2_mex,likec2_mex,err_mex,fevals_mex] = VestBMS_c2corrpostandlikec2adapt_mex(priorpdf2d,like_vis,like_vest); t_mex = toc;

% Correctness check
postright_c2_err = sum(abs(postright_c2(:)-postright_c2_mex(:)));
fprintf('Total error (postright_c2): %g\n', postright_c2_err);
if postright_c2_err > TolErr*numel(postright_c2);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in postright_c2.');
end
likec2_err = sum(abs(likec2(:)-likec2_mex(:)));
fprintf('Total error (likec2): %g\n', likec2_err);
if likec2_err > TolErr*numel(likec2);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in likec2.');
end
err_err = sum(abs(err(:)-err_mex(:)));
fprintf('Total error (err): %g\n', err_err);
if err_err > TolErr*numel(err);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in err.');
end
fevals_err = sum(abs(fevals(:)-fevals_mex(:)));
fprintf('Total error (fevals): %g\n', fevals_err);
if fevals_err > TolErr*numel(fevals);
	error('mexxer:tooLargeError','Correctness check failed. Error too large in fevals.');
end

% Runtime analysis
fprintf('Time for MATLAB code: %.3f s\n', t);
fprintf('Time for MEX file: %.3f s\n', t_mex);
fprintf('Speed gain: %.2f\n', t/t_mex);
