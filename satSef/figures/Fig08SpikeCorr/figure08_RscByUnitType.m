function [] = figure08_RscByUnitType(rscData, useNeuronTypes)
%FIGURE08_RSCBYUNITTYPE Summary of this function goes here
monkey = 'Da_Eu';
groupCols = {'condition','satCondition','outcome'};
rhocols = {'absRho','absRhoEst40','absRhoEst80'};

oPdfFile = 'figure08_RscByUnitType.pdf';
oExcelFile = 'figure08_RscByUnitType.xlsx';

figure
for un = 1:numel(useNeuronTypes)
    useNeuronType = useNeuronTypes{un};    
    if strcmp(useNeuronType,'ALL_NEURONS')
        idxRos = true(size(rscData,1),1);
        titleStr = 'All neurons';
    elseif strcmp(useNeuronType,'ERROR_NEURONS')
        idxRos = rscData.isSefErrorUnit == 1;
        titleStr = 'Neurons signalling error';
    elseif strcmp(useNeuronType,'OTHER_NEURONS')
        idxRos = rscData.isSefErrorUnit == 0;
        titleStr = 'Other neurons';
    end    
    rscSatConditionStats = grpstats(rscData(idxRos,['satCondition', rhocols]),'satCondition',{'mean','std','sem'});
    rscSatConditionStats.Properties.RowNames = {};
    
    rscOutcomesStats = grpstats(rscData(idxRos,[groupCols rhocols]),groupCols,{'mean','std','sem'});
    rscOutcomesStats = sortrows(rscOutcomesStats,{'outcome','satCondition'});
    rscOutcomesStats.Properties.RowNames = {};
    % display 3 groups of 2 bars each
    % data = [1 2;3 4;5 6]
    % quick:
    accClr = [1 0.2 0.2];
    fasClr = [0.2 1.0 0.2];
    
    grpColors = {accClr;fasClr};
    
    %         outcomes = {'Correct','ErrorChoice','ErrorTiming'}';
    %         sat = {'Accurate','Fast'};
    idxAccu = ismember(rscOutcomesStats.satCondition,'Accurate');
    idxFast = ismember(rscOutcomesStats.satCondition,'Fast');
    
    accuTbl = rscOutcomesStats(idxAccu,{'outcome','mean_absRho','sem_absRho'});
    fastTbl = rscOutcomesStats(idxFast,{'outcome','mean_absRho','sem_absRho'});
    
    
    subplot(1,numel(useNeuronTypes),un);
    
    [barCentersTbl, ~] = plotGroupBarsWithErrors(accuTbl.outcome,...
        [accuTbl.mean_absRho fastTbl.mean_absRho],...
        [accuTbl.sem_absRho fastTbl.sem_absRho],...
        grpColors);
    
    % Add boxes for Confidence interval
    % Acurate_Error_Timing ci/percentile 10/90 for 40 subsamples
    barCenter = barCentersTbl.ErrorTiming(1);
    idx = ismember(rscData.condition,'AccurateErrorTiming');
    ci = getCi(abs(rscData.rhoEst40(idx)));
    overplotBox(barCenter,ci,'k','-');
    
    % Fast_Error_Choice ci/percentile 10/90 for 80 subsamples
    idx = ismember(rscData.condition,'FastErrorChoice');
    barCenter = barCentersTbl.ErrorChoice(2);
    ci = getCi(abs(rscData.rhoEst80(idx)));
    overplotBox(barCenter,ci,'k',':');
    
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
    title(titleStr,'FontWeight','bold','Interpreter','none');
    drawnow
    
    %% DO 2-way ANOVA
    anovRes = satAnova(rscData(:,{'absRho','satCondition','outcome'}),'model','interaction',...
        'doMultiCompare',true,...
        'multiCompareDisplay','off',...
        'alpha',0.05);
    anovaTbl = anovRes.anovaTbl;
    satConditionByOutcome = sortrows(anovRes.satCondition_outcome(:,[1:6 9]),...
        {'levelName1','levelName2'},{'ascend','descend'});
    outcomeStats = rscOutcomesStats(:,[1,5,6,8,9,11,12]);
    satStats = rscSatConditionStats(:,[1,3,4,6,7,9,10]);
    satStats.Properties.VariableNames = outcomeStats.Properties.VariableNames;
    allStats = [satStats;outcomeStats];
    
    sheetName_prefix = [monkey useNeuronType];
    writetable(anovaTbl,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_anova']);
    writetable(satConditionByOutcome,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Interaction']);
    writetable(allStats,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Stats']);
        
end
ppretty([8,5],'XMinorTick','off');
gcf;
print(oPdfFile,'-dpdf')

%% Do 2-way anova and write excel file with multiple sheets
% do anova

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
