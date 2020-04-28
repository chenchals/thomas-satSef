function [] = figure_08()
%FIGURE_08: Generate figure 8 for satSef paper
%  
% see also: CREATESPIKECORRWITHSUBSAMPLING,

%% cretae datafile -- already done by running. Takes time
% CREATESPIKECORRWITHSUBSAMPLING
% [spkCorr] = createSpikeCorrWithSubSampling();
% creates the datafile.
% datafile : 
%% Load spkCorr data created in the above step
spkCorr = load('dataProcessed/satSefPaper/rscSubSampl1K_PostSaccade.mat');
spkCorr = spkCorr.spkCorr;

colNames = getColNamesToUse();
spkCorr = spkCorr(:,colNames);
spkCorr.Properties.VariableNames = regexprep(colNames,'_150ms','');

spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.epoch = spkCorr.alignedName;

%%
epoch = 'PostSaccade';
spkCorr = spkCorr(ismember(spkCorr.epoch,epoch),:);
rscTable = table();
rscTable.PairUid = spkCorr.Pair_UID;
rscTable.monkey = spkCorr.X_monkey;
rscTable.X_area = spkCorr.X_area;
rscTable.Y_area = spkCorr.Y_area;
rscTable.condition = spkCorr.condition;
rscTable.satCondition = spkCorr.satCondition;
rscTable.outcome = spkCorr.outcome;
rscTable.epoch = spkCorr.epoch;
rscTable.rho = spkCorr.rhoRaw;
rscTable.rhoEst40 = spkCorr.rhoEstRaw_nTrials_40;
rscTable.rhoEst80 = spkCorr.rhoEstRaw_nTrials_80;
% add isSefErrorNeuron
rscTable.isSefErrorUnit = abs(spkCorr.X_errGrade) > 1 | abs(spkCorr.X_rewGrade) > 1;
warning('off')
%%
outSpkCorr = table();
outSpkCorr.PairUid = rscTable.PairUid;
outSpkCorr.monkey = rscTable.monkey;
outSpkCorr.unitArea1 = rscTable.X_area;
outSpkCorr.unitArea2 = rscTable.Y_area;
outSpkCorr.satCondition = regexprep(rscTable.condition,{'Correct','Error.*'},'');
outSpkCorr.outcome = regexprep(rscTable.condition,{'Accurate','Fast'},'');
outSpkCorr.satOutcome = rscTable.condition;
outSpkCorr.rscObserved = rscTable.rho;
outSpkCorr.rscEstimated_40RandomTrials = rscTable.rhoEst40;
outSpkCorr.rscEstimated_80RandomTrials = rscTable.rhoEst80;
outSpkCorr = sortrows(outSpkCorr,{'unitArea1','unitArea2','outcome','satCondition','rscObserved'});
oExcelFile = 'fig08_data.xlsx';
writetable(outSpkCorr,oExcelFile,'UseExcel',true,'Sheet','Rsc_PostSaccade');

%% RSC by Unit Type
monkeys = {'Da_Eu'};
errorTypes = {{'ALL_NEURONS','ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08_RscByUnitType.pdf';
fig08RscMonkUnitType(rscTable,pdfFilename,monkeys,errorTypes,pdfFilename);

%%  by monkey...All neurons
monkeys = {'Da','Eu'};
monkUnitTypes = {{'ALL_NEURONS'}
                 {'ALL_NEURONS'}};
pdfFilename = 'fig08Suppl_RscMonkByAllUnits.pdf';    
fig08RscMonkUnitType(rscTable(rscTable.rho < 0,:),monkeys,monkUnitTypes,pdfFilename);

%%  by monkey...error and other neurons
monkeys = {'Da','Eu'};
monkUnitTypes = {{'ERROR_NEURONS','OTHER_NEURONS'}
                 {'ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08Suppl_RscMonkByUnitType.pdf';             
fig08RscMonkUnitType(rscTable,monkeys,monkUnitTypes,pdfFilename);

%% by monkey...FEF and SC pairs with SEF
%  Da-SEF-FEF, Da-SEF-SC, Eu-SEF-SC, Da_Eu-SEF-SC
monkeys = {'Da','Eu','Da_Eu'};
monkUnitTypes = {{'FEF','SC'}
                 {'SC'}
                 {'SC'}};
pdfFilename = 'fig08Suppl_RscByArea.pdf';
fig08RscMonkUnitType(rscTable,monkeys,monkUnitTypes,pdfFilename);

%% RSC by monk by area by errorType
%  Da-SEF-FEF-Error, Da-SEF-FEF-Other
%  Da-SEF-SC-Error, Da-SEF-SC-Other

monkeys = {'Da'};
unitAreas = {{'FEF','SC'}};
errorTypes = {{'ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08_RscByUnitType.pdf';
fig08RscMonkUnitType(rscTable,pdfFilename,'monkeys',monkeys,'unitAreas',unitAreas,'errorTypes',errorTypes);

end
% 
%                 Da-SEF-FEF, Da-SEF-SC, Eu-SEF-SC, Da_Eu-SEF-SC
%                Da-SEF-FEF-Error, Da-SEF-FEF-Other
%                Da-SEF-SC-Error, Da-SEF-SC-Other
%                Eu-SEF-SC-Error, Eu-SEF-SC-Other
%                Da_Eu-SEF-SC-Error, Da_Eu-SEF-SC-Other


%%
%%
function [colNames] = getColNamesToUse()
colNames = {
    'Pair_UID'
    'X_monkey'
    'X_unit'
    'Y_unit'
    'X_area'
    'Y_area'
    'X_errGrade'
    'X_rewGrade'
    'xSpkCount_win_150ms'
    'ySpkCount_win_150ms'
    'xMeanFr_spkPerSec_win_150ms'
    'yMeanFr_spkPerSec_win_150ms'
    'condition'
    'alignedName'
    'nTrials'
    'rhoRaw'
    'signifRaw_05'
    'rhoEstRaw_nTrials_40'
    'rhoEstSem_nTrials_40'
    'ci95_nTrials_40'
    'rhoRawInCi95_nTrials_40'
    'rhoEstRaw_nTrials_80'
    'rhoEstSem_nTrials_80'
    'ci95_nTrials_80'
    'rhoRawInCi95_nTrials_80'
    };
end