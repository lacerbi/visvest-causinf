function VestBMS_plotPSEdemo(data,bin,xb)
%VESTBMS_PLOTPSEDEMO

if nargin < 2 || isempty(bin); bin = 10; end
if nargin < 3; xb = []; end

xbin = [-37.5,-25,-17.5,-10,-3.75,0,3.75,10,17.5,25,37.5];
if isempty(xb); xb = xbin(bin); end

fontsize = 16;
plots = VestBMS_defaults('plots');  % Get default plot info

x = linspace(-25,25,512);

hold on;

mu_all = data.psyright_mu(:,bin);

% plotpsy(D,[1,6]);
plotpsy(data,[1,bin]);
plotpsy(data,[2,bin]);
plotpsy(data,[3,bin]);

set(gca,'TickDir','out');
box off;
xlabel('$s_\mathrm{vest}$','Interpreter','LaTeX','FontSize',fontsize);
ylabel('Pr(right)','Interpreter','TeX','FontSize',fontsize);
axis([x(1),x(end),0,1]);

set(gca,'XTick',[-30,-15,0,15,30],'XTickLabel',{'-30°','-15°','0°','15°','30°'},'YTick',0:0.5:1);

if max(mu_all) < 5
    text(5,0.5,'PSE','FontSize',fontsize);
    text(10,0.8,['$s_\mathrm{vis} = ' num2str(xb) '^{\circ}$'],'Interpreter','LaTeX','FontSize',fontsize);
else
    text(-21,0.5,'L/R PSE','FontSize',fontsize);
    text(-23.75,0.8,['$s_\mathrm{vis} = ' num2str(xb) '^{\circ}$'],'Interpreter','LaTeX','FontSize',fontsize);    
end

    function plotpsy(D,idx)
        
        iNoise = idx(1);
        plots.NoiseColors(iNoise,:);
        
        mu = D.psyright_mu(iNoise,idx(2));
        sigma = D.psyright_sigma(iNoise,idx(2));        
        y = 0.5 + 0.5*erf((x - mu)/(sqrt(2)*sigma));
        plot([mu mu],[0 1],'--k','LineWidth',1);
        plot(x,y,'-k','LineWidth',3,'Color',plots.NoiseColors(iNoise,:));
        scatter(mu,0.5,6^2,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
        
        
        

    end

end

