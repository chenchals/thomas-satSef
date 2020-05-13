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
rscTable.pval = spkCorr.pvalRaw;
rscTable.signif_05 = spkCorr.signifRaw_05;
rscTable.rhoEst40 = spkCorr.rhoEstRaw_nTrials_40;
rscTable.rhoEst80 = spkCorr.rhoEstRaw_nTrials_80;
% add isSefErrorNeuron
rscTable.isSefErrorUnit = abs(spkCorr.X_errGrade) > 1 | abs(spkCorr.X_rewGrade) > 1;
rscTable.sefVisGrade = spkCorr.X_visGrade;
rscTable.sefMoveGrade = spkCorr.X_moveGrade;
rscTable.isSefUnitVis = spkCorr.X_visGrade > 1 & spkCorr.X_moveGrade == 0;
rscTable.isSefUnitMove = spkCorr.X_visGrade == 0 & spkCorr.X_moveGrade > 1;
rscTable.isSefUnitVisMove = spkCorr.X_visGrade > 1 & spkCorr.X_moveGrade > 1;
rscTable.isSefUnitOther = spkCorr.X_visGrade <= 1 | spkCorr.X_moveGrade <= 1;
warning('off')
%%
outSpkCorr = table();
outSpkCorr.PairUid = rscTable.PairUid;
outSpkCorr.monkey = rscTable.monkey;
outSpkCorr.unitArea1 = rscTable.X_area;
outSpkCorr.unitArea2 = rscTable.Y_area;
outSpkCorr.sefVisGrade = rscTable.sefVisGrade;
outSpkCorr.sefMoveGrade = rscTable.sefMoveGrade;
outSpkCorr.isSefErrorUnit = rscTable.isSefErrorUnit;
outSpkCorr.isSefUnitVis = rscTable.isSefUnitVis;
outSpkCorr.isSefUnitMove = rscTable.isSefUnitMove;
outSpkCorr.isSefUnitVisMove = rscTable.isSefUnitVisMove;
outSpkCorr.isSefUnitOther = rscTable.isSefUnitOther;
outSpkCorr.satCondition = regexprep(rscTable.condition,{'Correct','Error.*'},'');
outSpkCorr.outcome = regexprep(rscTable.condition,{'Accurate','Fast'},'');
outSpkCorr.satOutcome = rscTable.condition;
outSpkCorr.rscObserved = rscTable.rho;
outSpkCorr.pvalObserved = rscTable.pval;
outSpkCorr.signif05 = rscTable.signif_05;
outSpkCorr.rscEstimated_40RandomTrials = rscTable.rhoEst40;
outSpkCorr.rscEstimated_80RandomTrials = rscTable.rhoEst80;
outSpkCorr = sortrows(outSpkCorr,{'unitArea1','unitArea2','outcome','satCondition','rscObserved'});
oExcelFile = 'fig08_data.xlsx';
writetable(outSpkCorr,oExcelFile,'UseExcel',true,'Sheet','Rsc_PostSaccade');

%% 04/29/2020
% JS ? We must have a look at signed Rsc values.
% All signed Rsc values for monkeys combined and separate, for error, and non-error neurons
% Only Rsc values > 0 for monkeys combined and separate, for error, and non-error neurons
% Only Rsc values < 0 for monkeys combined and separate, for error, and non-error neurons
% 
%% [Absolute|Signed|Positive|Negative] Rsc bar plots for error/non-error by monks
useMonkeys = {'Da_Eu','Da','Eu'};
useErrorTypes = {{'ALL_NEURONS','ERROR_NEURONS','OTHER_NEURONS'}};
rhoTypes = {'Absolute','Signed','Positive','Negative'};
for rt = 1:numel(rhoTypes)
    rhoType = rhoTypes{rt};
    evalin('base','addCiBox = false;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useErrorTypes)
    evalin('base','addCiBox = true;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useErrorTypes)
end

%% [Absolute|Signed|Positive|Negative] Rsc bar plots for FEF/SC pairs by monks
useMonkeys = {'Da_Eu','Da','Eu'};
useAreaTypes = {{'ALL_NEURONS','FEF','SC'}};
rhoTypes = {'Absolute','Signed','Positive','Negative'};
for rt = 1:numel(rhoTypes)
    rhoType = rhoTypes{rt};
    evalin('base','addCiBox = false;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useAreaTypes)
    evalin('base','addCiBox = true;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useAreaTypes)
end

%% [Absolute|Signed|Positive|Negative] Rsc bar plots for Vis/Mov/VisMove/Other by monks
useMonkeys = {'Da_Eu','Da','Eu'};
useFuncTypes = {{'ALL_NEURONS','VIS','MOV','VISMOV','OTHER'}};
rhoTypes = {'Absolute','Signed','Positive','Negative'};
for rt = 1:numel(rhoTypes)
    rhoType = rhoTypes{rt};
    evalin('base','addCiBox = false;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useFuncTypes)
    evalin('base','addCiBox = true;')
    doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useFuncTypes)
end

end
% 
%                 Da-SEF-FEF, Da-SEF-SC, Eu-SEF-SC, Da_Eu-SEF-SC
%                Da-SEF-FEF-Error, Da-SEF-FEF-Other
%                Da-SEF-SC-Error, Da-SEF-SC-Other
%                Eu-SEF-SC-Error, Eu-SEF-SC-Other
%                Da_Eu-SEF-SC-Error, Da_Eu-SEF-SC-Other


%%
function [] = doBarplotAndAnovaFor(rscTable,rhoType,useMonkeys,useUnitTypes)
% rhoType: [Absolute|Signed|Positive|Negative] 
if sum(contains(useUnitTypes{:},'ERROR')) > 0 
    midfix = 'ErrNoErr_'; 
elseif sum(strcmpi(useUnitTypes{:},'FEF')) > 0 || sum(strcmpi(useUnitTypes{:},'SC')) > 0
    midfix = 'FefSc_';
elseif sum(strcmpi(useUnitTypes{:},'VIS')) > 0 || sum(strcmpi(useUnitTypes{:},'MOV')) > 0
    midfix = 'VisMov_';
end

switch rhoType
    case 'Absolute'
        rscTable.rho = abs(rscTable.rho);
        rscTable.rhoEst40 = abs(rscTable.rhoEst40);
        rscTable.rhoEst80 = abs(rscTable.rhoEst80);
        basePdfFile = ['absoluteRsc_' midfix];
        titlePrefix = 'Absolute Rsc - ';
    case 'Signed'
        % No changes use Rho values as is
        basePdfFile = ['signedRsc_' midfix];
        titlePrefix = 'Signed Rsc - ';
    case 'Positive'
        rscTable = rscTable(rscTable.rho > 0,:);
        basePdfFile = ['positiveRsc_' midfix];
        titlePrefix = 'Positive Rsc - ';
    case 'Negative'        
        rscTable = rscTable(rscTable.rho < 0,:);
        basePdfFile = ['negativeRsc_' midfix];
        titlePrefix = 'Negative Rsc - ';
end

for m = 1:numel(useMonkeys)
    monkeys = useMonkeys(m);
    unitTypes = useUnitTypes;
    pdfFilename = [basePdfFile monkeys{1} '.pdf'];
    fig08RscMonkUnitType(rscTable,monkeys,unitTypes,pdfFilename,[titlePrefix monkeys{1}]);
    close all
end

end

%%
function [colNames] = getColNamesToUse()
colNames = {
    'Pair_UID'
    'X_monkey'
    'X_unit'
    'Y_unit'
    'X_area'
    'Y_area'
    'X_visGrade'
    'X_moveGrade'
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
    'pvalRaw'
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