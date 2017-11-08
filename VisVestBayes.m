function VisVestBayes()

Ngrid = 4001;
srange = linspace(-180,180,Ngrid)';
ds = srange(2)-srange(1);
xrange = srange';

% sigmazero_vest = 12; w_vest = 0.02;
sigmazero_vest = 5; w_vest = 0.08;

sigma_vest = VestBMS_sensoryNoise('D',srange,sigmazero_vest,w_vest);
plot(srange,sigma_vest); hold on;

% Wrapped normal likelihood
like = bsxfun_normpdf(xrange,srange,sigma_vest) + bsxfun_normpdf(xrange,srange+360,sigma_vest) + bsxfun_normpdf(xrange,srange-360,sigma_vest);
like = bsxfun(@rdivide,like,qtrapz(like,1)*ds);

[~,idx_center] = max(like,[],1);   % Center each likelihood around peak

for i = 1:Ngrid
    idx_shift = 0.5*(Ngrid-1)+idx_center(i)-1;
    svec = circshift(srange,idx_shift);
    bias(i) = qtrapz(like(:,i).*svec,1)*ds;
    %shat2(i) = qtrapz(like.*(svec.^2),1)*ds;
    %shatstd(i) = sqrt(shat2(i) - shat(i)^2);
end

plot(xrange,bias);













end