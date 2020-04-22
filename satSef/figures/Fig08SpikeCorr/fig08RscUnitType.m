function [anovaResultTbl] = fig08RscUnitType(rscData, useNeuronTypes)
%FIG08RSCUNITTYPE Summary of this function goes here
monkey = 'Da_Eu';
figDir = 'fig08';
if ~exist(figDir,'dir')
    mkdir(figDir);
end
oPdfFile = fullfile(figDir,'fig08_RscByUnitType.pdf');
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
for un = 1:numel(useNeuronTypes)
    useNeuronType = useNeuronTypes{un};
    subplot(1,numel(useNeuronTypes),un);
    anovaResultTbl.(useNeuronType) = doRscBarPlots(rscData,monkey,useNeuronType,oExcelFile);
end
ppretty([8,5],'XMinorTick','off');
gcf;
print(oPdfFile,'-dpdf')
save(oMatFile,'-v7.3','-struct','anovaResultTbl');

end


