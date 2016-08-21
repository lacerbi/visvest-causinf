% function VestBMS_c1postandlikec1qtrapz_test(N,K)

N = 401;
K = 301; 
% if nargin < 1 || isempty(N); N = 401; end
% if nargin < 2 || isempty(K); K = 301; end

x1 = randn(N,1,K);
x2 = randn(N,K);

tic; [a1,a2] = VestBMS_c1postandlikec1qtrapz_mat(x1,x2); t1 = toc;
tic; [b1,b2] = VestBMS_c1postandlikec1qtrapz(x1,x2); t2 = toc;

fprintf('=============================\n');
fprintf('Testing VestBMS_likec1qtrapz:\n');
fprintf('=============================\n');

fprintf('Time for MATLAB code: %.3f s\n', t1);
fprintf('Time for MEX file: %.3f s\n', t2);
fprintf('Speed gain: %.2f\n', t1/t2);
fprintf('Total error (1st output): %g\n', sum(abs(a1(:)-b1(:))));
fprintf('Total error (2nd output): %g\n', sum(abs(a2(:)-b2(:))));