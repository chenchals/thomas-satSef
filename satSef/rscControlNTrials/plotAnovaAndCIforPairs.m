%%
spkCorr = load('newRscWithConstantSubSampling.mat');
spkCorr = spkCorr.spkCorr;

colNames = getColNamesToUse();
spkCorr = spkCorr(:,colNames);
spkCorr.Properties.VariableNames = regexprep(colNames,'_150ms','');

conditions = {'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
              'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
    };

rscData = table();
rscData.condition = spkCorr.condition;
rscData.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
rscData.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
rscData.rho = spkCorr.rhoRaw;
rscData.rhoEst40 = spkCorr.rhoEstRaw_nTrials_40;
rscData.rhoEst80 = spkCorr.rhoEstRaw_nTrials_80;
rscData.absRho = abs(spkCorr.rhoRaw);
rscData.absRhoEst40 = abs(spkCorr.rhoEstRaw_nTrials_40);
rscData.absRhoEst80 = abs(spkCorr.rhoEstRaw_nTrials_80);

%% Get statistics for the specific rho values
groupCols = {'condition','satCondition','outcome'};
rhocols = {'absRho','absRhoEst40','absRhoEst80'};

rscSatConditionStats = grpstats(rscData(:,['satCondition', rhocols]),'satCondition',{'mean','std','sem'});
rscSatConditionStats.Properties.RowNames = {};

rscOutcomesStats = grpstats(rscData(:,[groupCols rhocols]),groupCols,{'mean','std','sem'});
rscOutcomesStats = sortrows(rscOutcomesStats,{'outcome','satCondition'});
rscOutcomesStats.Properties.RowNames = {};
% display 3 groups of 2 bars each
% data = [1 2;3 4;5 6]
% quick:
accClr = [1 0.2 0.2];
fasClr = [0.2 1.0 0.2];

grpColors = {accClr;fasClr};

outcomes = {'Correct','ErrorChoice','ErrorTiming'}';
sat = {'Accurate','Fast'};
idx = ismember(rscOutcomesStats.satCondition,'Accurate');
accuTbl = rscOutcomesStats(idx,{'outcome','mean_absRho','sem_absRho'});
idx = ismember(rscOutcomesStats.satCondition,'Fast');
fastTbl = rscOutcomesStats(idx,{'outcome','mean_absRho','sem_absRho'});
figure
[barCentersTbl, errBarHandles] = plotGroupedBarsWithErrs(accuTbl.outcome,...
    [accuTbl.mean_absRho fastTbl.mean_absRho],...
    [accuTbl.sem_absRho fastTbl.sem_absRho],...
    grpColors);

% Add boxes for Confidence interval
% Acurate_Error_Timing ci/percentile 10/90 for 40 subsamples
barCenter = barCentersTbl.ErrorTiming(1);
idx = ismember(rscData.condition,'AccurateErrorTiming');
ci = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci,'k','-');
%ptile = prctile(abs(rscData.rhoEst40(idx)),[10 90]);
%overplotBox(barCenter,ptile,'k','-');

% Fast_Error_Choice ci/percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastErrorChoice');
barCenter = barCentersTbl.ErrorChoice(2);
ci = getCi(abs(rscData.rhoEst80(idx)));
%overplotBox(barCenter,ci,'k',':');
%ptile = prctile(abs(rscData.rhoEst80(idx)),[10 90]);
overplotBox(barCenter,ptile,'k',':');

% Accurate_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'AccurateCorrect');
barCenter = barCentersTbl.Correct(1);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k','-');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci80,'k',':');

% Fast_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastCorrect');
barCenter = barCentersTbl.Correct(2);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k','-');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci80,'k',':');
ppretty([3.5,3])
%% DO 2-Factor ANOVA satCondition by outcome
temp = rscData;
temp.epoch = repmat({'PostSaccade'},size(rscData,1),1);
[anovRes] = doAnovaConditionByOutcome(temp,'absRho');

anovaTbl = anovRes.PostSaccade.anovaTbl
%satConditionTbl = anovRes.PostSaccade.satCondition(:,[1:6,9])
%outcomeTbl = anovRes.PostSaccade.outcome(:,[1:6,9])
satConditionByOutcome = sortrows(anovRes.PostSaccade.satCondition_outcome(:,[1:6 9]),{'levelName1','levelName2'},{'ascend','descend'})

outcomeStats = rscOutcomesStats(:,[1,5,6,8,9,11,12])
satStats = rscSatConditionStats(:,[1,3,4,6,7,9,10])
satStats.Properties.VariableNames = outcomeStats.Properties.VariableNames;
allStats = [satStats;outcomeStats]



%% Do 2-way anova
function [anovaConditionByOutcome] = doAnovaConditionByOutcome(spkCorr,useRho)
% Do 2-Factor anova for the following -
%    Factor1: condition
%    Factor 2: outcome
% The data colums for input anovaTbl from spkCorr data above:
%   yVals     : rho or rhoZ
%   condition : Factor-1 for 2-Factor anova values
%                 Fast|Accurate
%   outcome   : Factor-2 for 2-Factor anova values:
%                 Correct|ErrorChoice|ErrorTiming

epochs = {'PostSaccade'};
anovaConditionByOutcome = struct();
    for ep = 1:numel(epochs)
        epoch = epochs{ep};
        idx = ismember(spkCorr.epoch,epoch);
        anovaTbl = table();
        anovaTbl.rho = spkCorr{idx,useRho};
        % 2-way anova satcondition by outcome
        % Factor-1
        anovaTbl.satCondition = spkCorr{idx,'satCondition'};
        % Factor-2
        anovaTbl.outcome = spkCorr{idx,'outcome'};
        % do anova
        anovaModel = 'interaction'; %  [linear|interaction]
        multiCompareFlag = 1; % [0|1]
        alpha = 0.05;
        fn = epoch;
        anovaConditionByOutcome.(fn) = satAnova(anovaTbl,anovaModel,multiCompareFlag,alpha);
    end
end


%%
function [pH] = overplotBox(barCenter,yBeginEnd,edgeColor,lineStyle)
xBeginEnd = barCenter + [-0.125 0.125];
xVec = [xBeginEnd fliplr(xBeginEnd)];
yVec = [yBeginEnd;yBeginEnd];
yVec = yVec(:)';
pH = patch(xVec,yVec,[1 1 1],'FaceAlpha',0.0);
set(pH,'EdgeColor',edgeColor)
set(pH,'LineStyle',lineStyle)
set(pH,'LineWidth',0.5)

end
function [normCi] = getCi(vec)
    vecMean = nanmean(vec);
    vecSem = nanstd(vec)/sqrt(numel(vec));
    % compute t-statistic for 0.025, 0.975
    ts = tinv([0.025,0.975],numel(vec)-1);
    normCi = vecMean + vecSem*ts;
end
%%
function [colNames] = getColNamesToUse()
colNames = {
    'xSpkCount_win_150ms'
    'ySpkCount_win_150ms'
    'xMeanFr_spkPerSec_win_150ms'
    'yMeanFr_spkPerSec_win_150ms'
    'condition'
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