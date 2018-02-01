function [kl,kappa,kappa_prime] = vonmisesvswrapnormal()

fontsize = 16;

kappa_min = 1/(2*pi)^2;
kappa_max = 10^2;

kappa = exp(linspace(log(kappa_min),log(kappa_max),1e2));

nMix = 3;   % Use three mixture components for wrapped normals
x = linspace(-pi,pi,1e4);

for iter = 1:numel(kappa)

    y1 = vonmises(x,kappa(iter));
    [kappa_prime(iter),kl(iter)] = fminbnd(@(kappa_) kldiv_vm2wn(kappa_,y1,x,nMix),0.01,100);
    
    y2 = wrapnormal(x,kappa_prime(iter),nMix);
    plot(x,y1,'k','LineWidth',4); hold on;
    plot(x,y2,'r:','LineWidth',2);
    
end

xlim([-pi,pi]);
set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;
xlabel('x','FontSize',fontsize);
ylabel('Probability density','FontSize',fontsize);

hl = legend('Von Mises','Wrapped normal');
set(hl,'Location','NorthWest','FontSize',fontsize,'Box','off');

end


function y = wrapnormal(x,kappa,nMix)

if nargin < 3 || isempty(nMix); nMix = 3; end

sigma = 1/sqrt(kappa);

y = zeros(size(x));

for k = 1:nMix
    y = y + normpdf(x + (k-(nMix+1)/2)*2*pi,0,sigma);
end

dx = x(2)-x(1);
y = y / (qtrapz(y)*dx);

end


function y = vonmises(x,kappa)

y = exp(kappa*cos(x));

dx = x(2)-x(1);
y = y / (qtrapz(y)*dx);

end

function kl = kldiv_vm2wn(kappa_prime,y1,x,nMix)

y2 = wrapnormal(x,kappa_prime,nMix);

z = y1 .* log(y1 ./ y2);
z(~isfinite(z)) = 0;

dx = x(2)-x(1);
kl = qtrapz(z)*dx;

end