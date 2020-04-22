function [anovaResultTbl] = fig08SupplRscMonkUnitType(rscData, useNeuronTypes)
%FIG08SUPPLRSCMONKUNITTYPE Summary of this function goes here
monkey = 'Da_Eu';
groupCols = {'condition','satCondition','outcome'};
rhocols = {'absRho','absRhoEst40','absRhoEst80'};

oPdfFile = 'figure08Suppl_RscByUnitType.pdf';
oExcelFile = 'figure08Suppl_RscByUnitType.xlsx';
oMatFile = 'figure08Suppl_RscByUnitType.mat';
oFiles = {oPdfFile,oExcelFile,oMatFile};
for oF = oFiles
    if exist(oF{1},'file')
        delete(oF{1});
    end
end

fx_ANOVA2_roi = @(levelName1,levelName2,oCome) find(...
    ismember(levelName1,['Accurate_' oCome])...
    & ismember(levelName2,['Fast_' oCome])...
    );
anovaResultTbl = struct();
anovaTxt = {};
%% Suppl. Figure A: all neurons by monkey





figure
monkeys = {'Da','Eu'};
for mon = 1:numel(monkeys)
    monkey = monkeys{mon};
    neuronTypes = {'ALL'};
end
end
%%
function [anovaResultsTbl] = doPlot(useNeuronType)
nPlots = 3;


for un = 1:numel(neuronTypes)
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
    signifStrs = satConditionByOutcome.signifStr(roIds);
    
    satConditionByOutcome = [satConditionByOutcome(roIds,:);satConditionByOutcome(roIdsOther,:)];
    outcomeStats = rscOutcomesStats(:,[1,4,5,6,7,8,9,11,12]);
    satStats = rscSatConditionStats(:,[1,2,3,4,5,6,7,9,10]);
    satStats.Properties.VariableNames = outcomeStats.Properties.VariableNames;
    allStats = [satStats;outcomeStats];
    % write output to excel file
    sheetName_prefix = [monkey useNeuronType];
    writetable(anovaTbl,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_anova']);
    writetable(satConditionByOutcome,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Interaction']);
    writetable(allStats,oExcelFile,'UseExcel',true,'Sheet',[sheetName_prefix '_Stats']);
    
    %% Plot results/means
    % display 3 groups of 2 bars each
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
    
    outcomeLabels = accuTbl.outcome;
    roIds = cell2mat(cellfun(@(o) fx_ANOVA2_roi(satConditionByOutcome.levelName1,satConditionByOutcome.levelName2,o),...
        outcomeLabels,'UniformOutput',false));
    signifStrs = satConditionByOutcome.signifStr(roIds);
    
    
    [barCentersTbl, ~] = plotGroupBarsWithErrors(outcomeLabels,...
        [accuTbl.mean_absRho fastTbl.mean_absRho],...
        [accuTbl.sem_absRho fastTbl.sem_absRho],...
        signifStrs,...
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
    
    
    %% Anova results to show?
    anovaResultTbl.(useNeuronType).satConditionByOutcome = satConditionByOutcome;
    anovaResultTbl.(useNeuronType).anovaTbl = anovaTbl;
    anovaResultTbl.(useNeuronType).allStats = allStats;
    %tempTxt = getAnovaText(satConditionByOutcome,anovaTbl,allStats);
    
    
end
ppretty([8,5],'XMinorTick','off');
gcf;
print(oPdfFile,'-dpdf')
% save mat file of anova results
save(oMatFile,'-v7.3','-struct','anovaResultTbl');




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
