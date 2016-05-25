function pp = VestBMS_performanceCheck(data,subjs,pp)
%VESTBMS_PERFORMANCECHECK Measure subjects' performance in time

if nargin < 2 || isempty(subjs); subjs = {[1:11],[12:14]}; end
if nargin < 3; pp = []; end

if ~iscell(subjs); subjs = {subjs}; end

Npanels = numel(subjs);
wlen = 200;

for iPanel = 1:Npanels
    if Npanels > 1; subplot(1,Npanels,iPanel); end    
    Nsubjs = numel(subjs{iPanel});

    if isempty(pp)
        p = NaN(Nsubjs,1e4);

        for i = 1:Nsubjs
            nid = subjs{iPanel}(i)

            % Collapse performance across all conditions
            X = [];    
            % X = [X; data{nid}.X.bimodalall{3}];
            for j = 1:3; X = [X; data{nid}.X.bimodalall{j}]; end

            % Order by trial number
            [~,ord] = sort(X(:,1),'ascend');
            X = X(ord,:);

            % Remove ambigous trials
            X(X(:,3) == 0,:) = []; 

            tstart = X(1,1);
            tend = X(end,1);

            % Running window
            for j = tstart:tend-wlen
                win = any(bsxfun(@eq, X(:,1), j:j+wlen),2);
                Z = X(win,[3 5]);
                pcorrect = sum(Z(:,1).*Z(:,2) > 0)/size(Z,1);        
                p(i,j-tstart+1) = pcorrect;
            end
        end

        % Trim NaNs at the end
        n = find(any(isnan(p),1),1);
        p (:,n:end) = [];
    else
        p = pp{iPanel};
    end

    idx = 1:wlen:size(p,2);
    x = (1:size(idx,2))*wlen-wlen+1;
    y = nanmean(p(:,idx),1);
    yerr = nanstd(p(:,idx),1)./sqrt(sum(isfinite(p(:,idx)),1));
    fill([x fliplr(x)], [y + yerr, fliplr(y - yerr)],0.8*[1 1 1],'LineStyle','none'); hold on;
    plot(x,y,'k-','LineWidth',1);
    set(gca,'TickDir','out');
    set(gcf,'Color','w');
    box off;

    xlabel('Trials');
    ylabel('Fraction correct');
    
    if iPanel == 1
        title('Performance across trials (humans)');
    else
        title('Performance across trials (monkeys)');        
    end

    axis([1 max(x), 0.5 1]);
    
end

end