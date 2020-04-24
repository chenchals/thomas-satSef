function [anovaResultTbl] = doRscBarPlots(rscData,monkey,useNeuronType,oExcelFile)
anovaResultTbl = struct();
fx_ANOVA2_roi = @(levelName1,levelName2,oCome) find(...
    ismember(levelName1,['Accurate_' oCome])...
    & ismember(levelName2,['Fast_' oCome])...
    );

groupCols = {'condition','satCondition','outcome'};
rhocols = {'absRho','absRhoEst40','absRhoEst80'};

%useNeuronType = useNeuronTypes{un};
if strcmp(useNeuronType,'ALL_NEURONS')
    idxRos = true(size(rscData,1),1);
    titleStr = 'All neurons';
elseif strcmp(useNeuronType,'ERROR_NEURONS')
    idxRos = rscData.isSefErrorUnit == 1;
    titleStr = 'Error neurons';
elseif strcmp(useNeuronType,'OTHER_NEURONS')
    idxRos = rscData.isSefErrorUnit == 0;
    titleStr = 'Other neurons';
elseif strcmp(useNeuronType,'FEF')
    idxRos = ismember(rscData.Y_area, useNeuronType);
    titleStr = 'SEF-FEF';
elseif strcmp(useNeuronType,'SC')
    idxRos = ismember(rscData.Y_area, useNeuronType);
    titleStr = 'SEF-SC';
end

if sum(idxRos) == 0
    anovaResultTbl = [];
    return;
end

rscSatConditionStats = grpstats(rscData(idxRos,['satCondition', rhocols]),'satCondition',{'mean','std','sem'});
rscSatConditionStats.Properties.RowNames = {};

rscOutcomesStats = grpstats(rscData(idxRos,[groupCols rhocols]),groupCols,{'mean','std','sem'});
rscOutcomesStats = sortrows(rscOutcomesStats,{'outcome','satCondition'});
rscOutcomesStats.Properties.RowNames = {};

%% DO 2-way ANOVA
anovRes = satAnova(rscData(:,{'absRho','satCondition','outcome'}),'model','interaction',...
    'doMultiCompare',true,...
    'multiCompareDisplay','off',...
    'alpha',0.05);
anovaTbl = anovRes.anovaTbl;
satConditionByOutcome = sortrows(anovRes.satCondition_outcome(:,[1:6 9]),...
    {'levelName1','levelName2'},{'ascend','descend'});
% get row nos for comparisions of interest
roIds = cell2mat(cellfun(@(o) fx_ANOVA2_roi(satConditionByOutcome.levelName1,satConditionByOutcome.levelName2,o),...
    unique(rscData.outcome),'UniformOutput',false));
roIdsOther = setdiff([1:size(satConditionByOutcome,1)]',roIds);

satConditionByOutcome = [satConditionByOutcome(roIds,:);satConditionByOutcome(roIdsOther,:)];
outcomeStats = rscOutcomesStats(:,[1,4,5,6,7,8,9,11,12]);
satStats = rscSatConditionStats(:,[1,2,3,4,5,6,7,9,10]);
satStats.Properties.VariableNames = outcomeStats.Properties.VariableNames;
allStats = [satStats;outcomeStats];
% add CI table
temp = zeros(size(allStats,1),1);
allStats.loCI_40 = temp;
allStats.hiCI_40 = temp;
allStats.loCI_80 = temp;
allStats.hiCI_80 = temp;


%% Plot results/means
% display 3 groups of 2 bars each

idxAccu = ismember(rscOutcomesStats.satCondition,'Accurate');
idxFast = ismember(rscOutcomesStats.satCondition,'Fast');

accuTbl = rscOutcomesStats(idxAccu,{'outcome','mean_absRho','sem_absRho'});
fastTbl = rscOutcomesStats(idxFast,{'outcome','mean_absRho','sem_absRho'});

outcomeLabels = accuTbl.outcome;
roIds = cell2mat(cellfun(@(o) fx_ANOVA2_roi(satConditionByOutcome.levelName1,satConditionByOutcome.levelName2,o),...
    outcomeLabels,'UniformOutput',false));
signifStrs = satConditionByOutcome.signifStr(roIds);


[barCentersTbl, ~] = plotGroupBarsWithErrors(outcomeLabels,...
    [accuTbl.mean_absRho fastTbl.mean_absRho],...
    [accuTbl.sem_absRho fastTbl.sem_absRho],...
    signifStrs);

% Add boxes for Confidence interval
% Acurate_Error_Timing ci/percentile 10/90 for 40 subsamples
barCenter = barCentersTbl.ErrorTiming(1);
idx = ismember(rscData.condition,'AccurateErrorTiming');
ci = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci,'k','-');
idx = ismember(allStats.condition,'AccurateErrorTiming');
allStats.loCI_40(idx) = ci(1);
allStats.hiCI_40(idx) = ci(2);


% Fast_Error_Choice ci/percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastErrorChoice');
barCenter = barCentersTbl.ErrorChoice(2);
ci = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci,'k',':');
idx = ismember(allStats.condition,'FastErrorChoice');
allStats.loCI_80(idx) = ci(1);
allStats.hiCI_80(idx) = ci(2);

% Accurate_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'AccurateCorrect');
barCenter = barCentersTbl.Correct(1);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k','-');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci80,'k',':');
idx = ismember(allStats.condition,'AccurateCorrect');
allStats.loCI_40(idx) = ci40(1);
allStats.hiCI_40(idx) = ci40(2);
allStats.loCI_80(idx) = ci80(1);
allStats.hiCI_80(idx) = ci80(2);

% Fast_Correct percentile 10/90 for 80 subsamples
idx = ismember(rscData.condition,'FastCorrect');
barCenter = barCentersTbl.Correct(2);
ci40 = getCi(abs(rscData.rhoEst40(idx)));
overplotBox(barCenter,ci40,'k','-');
ci80 = getCi(abs(rscData.rhoEst80(idx)));
overplotBox(barCenter,ci80,'k',':');
idx = ismember(allStats.condition,'FastCorrect');
allStats.loCI_40(idx) = ci40(1);
allStats.hiCI_40(idx) = ci40(2);
allStats.loCI_80(idx) = ci80(1);
allStats.hiCI_80(idx) = ci80(2);

title([monkey '--' titleStr],'FontWeight','bold','Interpreter','none');
drawnow

% Anova results to show?
anovaResultTbl.satConditionByOutcome = satConditionByOutcome;
anovaResultTbl.anovaTbl = anovaTbl;
anovaResultTbl.allStats = allStats;

% write output to excel file
writeToExcel = false;
if writeToExcel
    sheetName_prefix = [monkey useNeuronType]; %#ok<UNRCH>
    writetable(anovaTbl,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_anova']);
    writetable(satConditionByOutcome,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Interaction']);
    writetable(allStats,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Stats']);
end

end

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


