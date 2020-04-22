function [anovaResultTbl] = fig08RscUnitType(rscData, useNeuronTypes)
%FIG08RSCUNITTYPE Summary of this function goes here
monkey = 'Da_Eu';

oPdfFile = 'figure08_RscByUnitType.pdf';
oExcelFile = 'figure08_RscByUnitType.xlsx';
oMatFile = 'figure08_RscByUnitType.mat';
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


