function [barCentersTbl, errBarHandles] = plotGroupedBarsWithErrs(cellArrCategories,yMean,yErr,varargin) %#ok<*INUSL>

   if numel(varargin) == 0
       accClr = [1 0.2 0.2];
       fasClr = [0.2 1.0 0.2];

       grpColors = {accClr;fasClr};
   else
       grpColors = varargin{1};
   end

    % find minRho for negative Rho 
    % manually was [-0.2 0.45]
    yLims = [0 0.18];
    
    temp = yMean(:);
    if max(temp) > max(yLims)
        yLims = [0 0.26]; % manually checked...
    end
    
    %x = categorical(xData(:,1));
    x = 1:numel(cellArrCategories);
    hBar = bar(x,yMean,'FaceAlpha',0.6,'BarWidth',0.9);
    if ~isempty(grpColors)
        set(hBar(1),'FaceColor',grpColors{1}); 
        set(hBar(2),'FaceColor',grpColors{2});
    end
    set(gca,'xticklabels',cellArrCategories,'TickLabelInterpreter','none','XTickLabelRotation',20);
    %set(gca,'YGrid','on')
    ylabel('Mean r_{sc} \pm SEM (abs.)');
    set(gca,'Ylim',yLims);
    set(gca,'FontWeight','bold','FontSize',8);
    hold on
    for k1 = 1:size(yMean,2)
        barCenters(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]'); %#ok<AGROW>
    end
    errBarHandles = errorbar(barCenters,yMean',yErr','Marker','none','MarkerSize',5,'LineStyle',...
        'none','LineWidth',0.5,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');

    
    barCentersTbl = array2table(barCenters,'VariableNames',cellArrCategories);
    barCentersTbl.Properties.RowNames = arrayfun(@(x) sprintf('Group#%d',x),1:size(barCentersTbl,1),'UniformOutput',false);
    
    
%     %%add text values to plot
%     if isPlusRhos
%         set(errBarHandles,'YPositiveDelta',[]);
%         y = 0.2;
%         va = 'bottom';
%     else
%         set(errBarHandles,'YNegativeDelta',[]);
%          y = -0.2;
%         va = 'top';
%    end    
%     
%     
%     [txtPairArea,txtAcc,txtFast] = arrayfun(@(x) deal(sprintf('%14s',xData{x,1}),...
%         sprintf('%0.2f\\pm%0.2f %4d',yMean(x,1),yErr(x,1),nPairs(x,1)),...
%         sprintf('%0.2f\\pm%0.2f %4d',yMean(x,2),yErr(x,2),nPairs(x,2))),...
%         (1:size(yMean,1))','UniformOutput',false);
%     txtPairArea = [sprintf('%14s','AREA');txtPairArea];
%     txtAcc = [sprintf('%14s','\mu\pmSEM   n');txtAcc];
%     txtFast = [sprintf('%14s','\mu\pmSEM   n');txtFast];
%     
% 
% 
%     h_pa = text(barCenters(2,2)+0.5,y,txtPairArea,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7);
%     ext = get(h_pa,'Extent');
%     pos = get(h_pa,'Position');
%     h_txt1 = text(pos(1)+ext(3),y,txtAcc,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7,'color',accColor);
%     ext = get(h_txt1,'Extent');
%     pos = get(h_txt1,'Position');
%     h_txt2 = text(pos(1)+ext(3)+0.1,y,txtFast,'HorizontalAlignment','center','VerticalAlignment',va,'FontSize',7,'color',fastColor); %#ok<NASGU>
    
end

