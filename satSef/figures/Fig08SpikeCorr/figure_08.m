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
spkCorr = load('dataProcessed/satSefPaper/newRscSubSampling1K_AllEpochs.mat');
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
rscPostSaccade = table();
rscPostSaccade.monkey = spkCorr.X_monkey;
rscPostSaccade.Y_area = spkCorr.Y_area;
rscPostSaccade.condition = spkCorr.condition;
rscPostSaccade.satCondition = spkCorr.satCondition;
rscPostSaccade.outcome = spkCorr.outcome;
rscPostSaccade.epoch = spkCorr.epoch;
rscPostSaccade.rho = spkCorr.rhoRaw;
rscPostSaccade.rhoEst40 = spkCorr.rhoEstRaw_nTrials_40;
rscPostSaccade.rhoEst80 = spkCorr.rhoEstRaw_nTrials_80;
rscPostSaccade.absRho = abs(spkCorr.rhoRaw);
rscPostSaccade.absRhoEst40 = abs(spkCorr.rhoEstRaw_nTrials_40);
rscPostSaccade.absRhoEst80 = abs(spkCorr.rhoEstRaw_nTrials_80);
% add isSefErrorNeuron
rscPostSaccade.isSefErrorUnit = abs(spkCorr.X_errGrade) > 1 | abs(spkCorr.X_rewGrade) > 1;
warning('off')

%% RSC by Unit Type
monkeys = {'Da_Eu'};
errorTypes = {{'ALL_NEURONS','ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08_RscByUnitType.pdf';
fig08RscMonkUnitType(rscPostSaccade,pdfFilename,'monkeys',monkeys,'errorTypes',errorTypes);

%%  by monkey...All neurons
monkeys = {'Da','Eu'};
monkUnitTypes = {{'ALL_NEURONS'}
                 {'ALL_NEURONS'}};
pdfFilename = 'fig08Suppl_RscMonkByAllUnits.pdf';             
fig08RscMonkUnitType(rscPostSaccade,monkeys,monkUnitTypes,pdfFilename);

%%  by monkey...error and other neurons
monkeys = {'Da','Eu'};
monkUnitTypes = {{'ERROR_NEURONS','OTHER_NEURONS'}
                 {'ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08Suppl_RscMonkByUnitType.pdf';             
fig08RscMonkUnitType(rscPostSaccade,monkeys,monkUnitTypes,pdfFilename);

%% by monkey...FEF and SC pairs with SEF
%  Da-SEF-FEF, Da-SEF-SC, Eu-SEF-SC, Da_Eu-SEF-SC
monkeys = {'Da','Eu','Da_Eu'};
monkUnitTypes = {{'FEF','SC'}
                 {'SC'}
                 {'SC'}};
pdfFilename = 'fig08Suppl_RscByArea.pdf';
fig08RscMonkUnitType(rscPostSaccade,monkeys,monkUnitTypes,pdfFilename);

%% RSC by monk by area by errorType
%  Da-SEF-FEF-Error, Da-SEF-FEF-Other
%  Da-SEF-SC-Error, Da-SEF-SC-Other

monkeys = {'Da'};
unitAreas = {{'FEF','SC'}};
errorTypes = {{'ERROR_NEURONS','OTHER_NEURONS'}};
pdfFilename = 'fig08_RscByUnitType.pdf';
fig08RscMonkUnitType(rscPostSaccade,pdfFilename,'monkeys',monkeys,'unitAreas',unitAreas,'errorTypes',errorTypes);

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
    'nSubSamples'
    'nTrials4SubSample'
    'rhoRaw'
    'signifRaw_05'
    'rhoEstRaw'
    'rhoEstSem'
    'prctile_10_90'
    'ci95'
    'rhoRawInCi95'
    'rhoRawInPrctile_10_90'
    'rhoEstRaw_nTrials_40'
    'rhoEstSem_nTrials_40'
    'prctile_10_90_nTrials_40'
    'ci95_nTrials_40'
    'rhoRawInCi95_nTrials_40'
    'rhoRawInPrctile_10_90_nTrials_40'
    'rhoEstRaw_nTrials_80'
    'rhoEstSem_nTrials_80'
    'prctile_10_90_nTrials_80'
    'ci95_nTrials_80'
    'rhoRawInCi95_nTrials_80'
    'rhoRawInPrctile_10_90_nTrials_80'
    };
end