function [anovaResultTbl] = fig08RscMonkUnitType(rscData, monkeys, cellArrUnitTypes,pdfFilename,titleStr)
%FIG08RSCMONKUNITTYPE Summary of this function goes here

% check if if you need to add CI box to the barplots
doCiBox = true;

if evalin('base','exist(''addCiBox'',''var'')')
    doCiBox = evalin('base','addCiBox');
end
if doCiBox
    figDir = '_fig08_With_CI_Box';
else
    figDir = '_fig08_No_CI_Box';

end

if ~exist(figDir,'dir')
    mkdir(figDir);
end

oPdfFile =  fullfile(figDir,pdfFilename);
oExcelFile = regexprep(oPdfFile,'pdf$','xlsx');
oMatFile = regexprep(oPdfFile,'pdf$','mat');
oFiles = {oPdfFile,oExcelFile,oMatFile};
for oF = oFiles
    if exist(oF{1},'file')
        delete(oF{1});
    end
end

anovaResultTbl = struct();
%figure
nSubplots = sum(cellfun(@(x) numel(x),cellArrUnitTypes));
plotNo = 0;
H_axes = figTemplate(1,nSubplots);
for mon = 1:numel(monkeys)
    monkey = monkeys{mon};
    if strcmp(monkey,'Da_Eu')
        monkIdx = true(size(rscData,1),1);
    else
        monkIdx = ismember(rscData.monkey,monkey(1));
    end
    unitTypes = cellArrUnitTypes{mon};
    for un = 1:numel(unitTypes)
        useNeuronType = unitTypes{un};
        monkData = rscData(monkIdx,:);
        plotNo = plotNo + 1;
        
        %subplot(1,nSubplots,plotNo);
        axes(H_axes(plotNo));
        
        temp = doRscBarPlots(monkData,monkey,useNeuronType,oExcelFile);
        if ~isempty(temp)
            anovaResultTbl.(monkey).(useNeuronType) = temp;
        else
            plotNo = plotNo - 1;
        end
    end
end
gcf;
ha = annotation('textbox','String',titleStr,...
    'Position',[0.25,0.95,0.50,0.05],...
    'HorizontalAlignment','left','LineStyle','none',...
    'Interpreter','none','FontWeight','bold','FontSize',12);
set(H_axes,'Box','off','TickDir','out','XMinorTick','off','YMinorTick','on');
set(gcf,'Position',[120 120 1400 900]);
set(gcf,'PaperOrientation','landscape')
ppretty([5 7],'XMinorTick','off')
drawnow
print(oPdfFile,'-dpdf')

%save(oMatFile,'-v7.3','-struct','anovaResultTbl');

end


function [hA] = figTemplate(ros,cols)
%%
figure;
set(gcf,'Color',[1 1 1]);
hA = tight_subplot(ros, cols, [.08 .08],[0.1 0.1],[0.1 0.1]);
for ii = 1:numel(hA)
    axes(hA(ii))
    title(sprintf('Plot %d',ii));
end
end