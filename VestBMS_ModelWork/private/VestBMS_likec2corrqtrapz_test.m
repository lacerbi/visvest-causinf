% function VestBMS_likec2corrqtrapz_test(N,K)

N = 201;
K = 251;

if isempty(N); N = 201; end
if isempty(K); K = 251; end

p = randn(N,N);
x1 = randn(N,K);
x2 = randn(N,1,K);

a = []; b = []; t1 = []; t2 = [];
tic; a = VestBMS_likec2corrqtrapz_mat(p,x1,x2); t1 = toc;
tic; b = VestBMS_likec2corrqtrapz(p,x1,x2); t2 = toc;

fprintf('=============================\n');
fprintf('Testing VestBMS_likec2corrqtrapz:\n');
fprintf('=============================\n');

if ~isempty(t1); fprintf('Time for MATLAB code: %.3f s\n', t1); end
if ~isempty(t2); fprintf('Time for MEX file: %.3f s\n', t2); end
if ~isempty(t1) && ~isempty(t2); fprintf('Speed gain: %.2f\n', t1/t2); end
if ~isempty(a) && ~isempty(b)
    fprintf('Total error: %g\n', sum(abs(a(:)-b(:))));
end