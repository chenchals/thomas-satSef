function [anovaResultTbl] = fig08SupplRscMonkUnitType(rscData, monkeys, useNeuronTypes)
%FIG08SUPPLRSCMONKUNITTYPE Summary of this function goes here
if numel(useNeuronTypes) == 1
    midfix = '_A_';
elseif numel(useNeuronTypes) == 2 && sum(contains(useNeuronTypes,'FEF'))
    midfix = '_A2_';
elseif numel(useNeuronTypes) == 2
    midfix = '_B_';
else
    midfix = '_';
end
figDir = 'fig08';
if ~exist(figDir,'dir')
    mkdir(figDir);
end

oPdfFile =  fullfile(figDir,['fig08Suppl' midfix 'RscByUnitType.pdf']);
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
nSubplots = numel(monkeys)*numel(useNeuronTypes);
plotNo = 0;
for mon = 1:numel(monkeys)
    monkey = monkeys{mon};
    for un = 1:numel(useNeuronTypes)
        useNeuronType = useNeuronTypes{un};
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