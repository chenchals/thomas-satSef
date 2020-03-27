%load r_sc data
% Use onkly SEF_FEF and SEF_SC pairs
load('dataProcessed/satSefPaper/analysis/spkCorr/summary/SAT_SEF_StaticRscAllPairs.mat', 'SEF_FEF','SEF_SC');
useCols = {'pairAreas','condition','alignedName','nTrials','rhoRaw_150ms'};
% useCols
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
% recode columns for anova
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.epoch = spkCorr.alignedName;
% take absolute value of correlation
spkCorr.rhoR = abs(spkCorr.rhoRaw_150ms);

%
maxTrials = max(spkCorr.nTrials);
minTrials = min(spkCorr.nTrials);

% conditions
satConditions = {'Fast','Accurate'};
epochs = {'PostSaccade'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};

