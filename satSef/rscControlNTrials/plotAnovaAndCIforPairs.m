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
rscDataStats = grpstats(rscData(:,[groupCols rhocols]),groupCols,{'mean','sem'});
rscDataStats = sortrows(rscDataStats,{'outcome','satCondition'});
rscDataStats.Properties.RowNames = {};
% display 3 groups of 2 bars each
% data = [1 2;3 4;5 6]
%rsc groups
% rscDataStats =
% 
%   6×10 table
% 
%           condition          satCondition       outcome       GroupCount    mean_absRho    sem_absRho    mean_absRhoEst40    sem_absRhoEst40    mean_absRhoEst80    sem_absRhoEst80
%     _____________________    ____________    _____________    __________    ___________    __________    ________________    _______________    ________________    _______________
% 
%     'AccurateCorrect'         'Accurate'     'Correct'           258         0.082981      0.0042224         0.081999           0.0042545           0.082866           0.0042311   
%     'FastCorrect'             'Fast'         'Correct'           258         0.082178      0.0046837         0.080563           0.0046447           0.081294           0.0047255   
%     'AccurateErrorChoice'     'Accurate'     'ErrorChoice'       258           0.1174      0.0054808          0.11567           0.0054272            0.11709           0.0054492   
%     'FastErrorChoice'         'Fast'         'ErrorChoice'       258          0.12106      0.0068834          0.12054           0.0069491            0.12085           0.0068788   
%     'AccurateErrorTiming'     'Accurate'     'ErrorTiming'       258          0.10216      0.0054097          0.10115           0.0053992            0.10155           0.0054007   
%     'FastErrorTiming'         'Fast'         'ErrorTiming'       258           0.1621      0.0082477          0.15845           0.0081643            0.15909           0.0081793   

% quick:
accClr = [1 0.2 0.2];
fasClr = [0.2 1.0 0.2];

grpColors = {accClr;fasClr};

outcomes = {'Correct','ErrorChoice','ErrorTiming'}';
sat = {'Accurate','Fast'};
idx = ismember(rscDataStats.satCondition,'Accurate');
accuTbl = rscDataStats(idx,{'outcome','mean_absRho','sem_absRho'});
idx = ismember(rscDataStats.satCondition,'Fast');
fastTbl = rscDataStats(idx,{'outcome','mean_absRho','sem_absRho'});
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
overplotBox(barCenter,ci,'k');
ptile = prctile(abs(rscData.rhoEst40(idx)),[10 90]);
%overplotBox(barCenter,ptile,'k');

% Fast_Error_Choice ci/percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastErrorChoice');
barCenter = barCentersTbl.ErrorChoice(2);
ci = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci,'k');
ptile = prctile(abs(rscData.rhoEst80(idx)),[10 90]);
%overplotBox(barCenter,ptile,'k');

% Accurate_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'AccurateCorrect');
barCenter = barCentersTbl.Correct(1);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
%overplotBox(barCenter,ci80,'k');

% Fast_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastCorrect');
barCenter = barCentersTbl.Correct(2);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci80,'k');




% xBeginEnd = barCentersTbl.ErrorTiming(1) + [-0.125 0.125];
% yBeginEnd = ci;
% xVec = [xBeginEnd fliplr(xBeginEnd)];
% yVec = [yBeginEnd;yBeginEnd];
% yVec = yVec(:)';
% pH = patch(xVec,yVec,[0.9 0.9 0.9],'FaceAlpha',0.1)

%%
function [pH] = overplotBox(barCenter,yBeginEnd,edgeColor)
xBeginEnd = barCenter + [-0.125 0.125];
xVec = [xBeginEnd fliplr(xBeginEnd)];
yVec = [yBeginEnd;yBeginEnd];
yVec = yVec(:)';
pH = patch(xVec,yVec,[1 1 1],'FaceAlpha',0.0);
set(pH,'EdgeColor',edgeColor)
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