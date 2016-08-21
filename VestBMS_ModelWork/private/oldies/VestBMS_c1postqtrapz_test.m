function VestBMS_c1postqtrapz_test(N,K)

if nargin < 1 || isempty(N); N = 401; end
if nargin < 2 || isempty(K); K = 301; end

x1 = randn(N,1,K);
x2 = randn(N,K);

tic; a = VestBMS_PostRight(bsxfun(@times,x1,x2)); t1 = toc;
tic; b = VestBMS_c1postqtrapz(x1,x2); t2 = toc;

fprintf('=============================\n');
fprintf('Testing VestBMS_c1postqtrapz:\n');
fprintf('=============================\n');

fprintf('Time for MATLAB code: %.3f s\n', t1);
fprintf('Time for MEX file: %.3f s\n', t2);
fprintf('Speed gain: %.2f\n', t1/t2);
fprintf('Total error: %g\n', sum(abs(a(:)-b(:))));

end