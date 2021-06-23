%load r_sc data
% Use onkly SEF_FEF and SEF_SC pairs
load('dataProcessed/satSefPaper/analysis/spkCorr/summary/SAT_SEF_StaticRscAllPairs.mat', 'SEF_FEF','SEF_SC');

useCols = {'pairAreas','condition','alignedName','nTrials','rhoRaw_150ms','signifRaw_05_150ms','signifRaw_01_150ms'};
% useCols
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
% recode columns for anova
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.epoch = spkCorr.alignedName;
% take absolute value of correlation
spkCorr.absRho = abs(spkCorr.rhoRaw_150ms);
spkCorr.rho_squared = spkCorr.rhoRaw_150ms.^2;
% signifAt05
spkCorr.signifAt05 = spkCorr.signifRaw_05_150ms;
% signifAt01
spkCorr.signifAt01 = spkCorr.signifRaw_01_150ms;
% conditions
satConditions = {'Fast','Accurate'};
epochs = {'PostSaccade'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
%
maxTrials = max(spkCorr.nTrials);
minTrials = min(spkCorr.nTrials);

%% JS counts of Rsc by condition - CS added different significance levels
rscCountsByCond = grpstats(spkCorr(:,{'outcome','satCondition','condition','nTrials'}),...
    {'outcome','satCondition','condition'},{'min','median','max','mode'});
rscCountsByCond.Properties.RowNames = {};
rscCountsByCond.countAll = rscCountsByCond.GroupCount;
rscCountsByCond.GroupCount = [];
% only signif05
temp =  grpstats(spkCorr(spkCorr.signifAt05,{'condition'}),'condition','count');
temp.Properties.RowNames = {};
rscCountsByCond = join(rscCountsByCond,temp);
rscCountsByCond.countAtSignif05 = rscCountsByCond.GroupCount;
rscCountsByCond.GroupCount = [];

% only signif01
temp =  grpstats(spkCorr(spkCorr.signifAt01,{'condition'}),'condition','count');
temp.Properties.RowNames = {};
rscCountsByCond = join(rscCountsByCond,temp);
rscCountsByCond.countAtSignif01 = rscCountsByCond.GroupCount;
rscCountsByCond.GroupCount = [];
%% rsc significant at p<0.05
rowIds = find(spkCorr.signifAt05 > 0 & ismember(spkCorr.epoch,'PostSaccade'));
rscNtrialsByCondSig05 =  getGroupStats(rowIds,spkCorr)


%% rsc non significant 
rowIds = find(spkCorr.signifAt05 == 0 & ismember(spkCorr.epoch,'PostSaccade'));
rscNtrialsByCondNonSig =  getGroupStats(rowIds,spkCorr)

%% rsc all
rowIds = find(spkCorr.signifAt05 >= 0 & ismember(spkCorr.epoch,'PostSaccade'));
rscNtrialsByCondAll =  getGroupStats(rowIds,spkCorr)




%% TR nTrials distributions by condition
bins = minTrials:maxTrials;
binEdges = minTrials-0.5:maxTrials+0.5;
fx_trialDist = @(x) histcounts(x,binEdges);
nTrlsDistrib


%%
function [statsTbl] = getGroupStats(rowIds,spkCorr)
    temp = grpstats(spkCorr(rowIds,{'outcome','satCondition','condition','nTrials'}),...
        {'outcome','satCondition','condition'},{'min','max','mode'});
    temp.Properties.RowNames = {};
    temp.countAll = temp.GroupCount;
    temp.GroupCount = [];
    temp1 = temp(:,3:end);
    temp = grpstats(spkCorr(rowIds,{'condition','absRho','rho_squared'}),...
        {'condition'},{'mean','std'});
    temp.Properties.RowNames = {};
    temp.GroupCount = [];
    statsTbl = temp1;
    statsTbl = join(statsTbl,temp);
end

