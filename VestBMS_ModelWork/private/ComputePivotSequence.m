function sbar = ComputePivotSequence(sigmazero,w,gamma)

sbar(1) = 0;

i = 1;
while sbar(i) < 45
    i = i + 1;
    sbar(i) = (sbar(i-1)*gamma + sqrt(sbar(i-1)^2*gamma^2 + (1-gamma)*(4*gamma*sigmazero^2/w^2 + (1+gamma)*sbar(i-1)^2)))/(1-gamma); 
end

srange = linspace(-45,45,2001);
hold off;
sbar = [-sbar(end:-1:2),sbar];
plot(sbar, sqrt(sigmazero^2 + w^2.*sbar.^2),'r','LineWidth',3);
hold on;
y = sqrt(sigmazero^2 + w^2.*srange.^2);
plot(srange, y,'k','LineWidth',1);
box off;
axis([-45,45,0,max(y)]);
set(gcf,'Color','w');
set(gca,'TickDir','out');

end