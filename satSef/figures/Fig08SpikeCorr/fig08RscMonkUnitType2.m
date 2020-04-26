function [anovaResultTbl] = fig08RscMonkUnitType(rscData, monkeys, cellArrUnitTypes,pdfFilename)
%FIG08RSCMONKUNITTYPE Summary of this function goes here

figDir = 'fig08_';
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
figure
nSubplots = sum(cellfun(@(x) numel(x),cellArrUnitTypes));
plotNo = 0;
for mon = 1:numel(monkeys)
    monkey = monkeys{mon};
    unitTypes = cellArrUnitTypes{mon};
    for un = 1:numel(unitTypes)
        useNeuronType = unitTypes{un};
        monkData = rscData(contains(rscData.monkey,monkey(1)),:);
        plotNo = plotNo + 1;
        subplot(1,nSubplots,plotNo);
        temp = doRscBarPlots(monkData,monkey,useNeuronType,oExcelFile);
        if ~isempty(temp)
        anovaResultTbl.(monkey).(useNeuronType) = temp;
        else
            plotNo = plotNo - 1;
        end
    end
end
ppretty([8,5],'XMinorTick','off');
gcf;
print(oPdfFile,'-dpdf')
save(oMatFile,'-v7.3','-struct','anovaResultTbl');



end