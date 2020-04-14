%%
spkCorr = load('newRscWithConstantSubSampling.mat');
spkCorr = spkCorr.spkCorr;

colNames = getColNamesToUse();
spkCorr = spkCorr(:,colNames);
spkCorr.Properties.VariableNames = regexprep(colNames,'_150ms','');

conditions = {'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
              'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
    };
spkCorr.rhoRawAbs = abs(spkCorr.rhoRaw);
spkCorr.rhoRawSqr = spkCorr.rhoRaw.^2;

spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
% add isSefErrorNeuron
spkCorr.isSefErrorUnit = abs(spkCorr.X_errGrade) > 1 & abs(spkCorr.X_rewGrade) > 1;

%% what is the distribution of positive and negative rhos overall
% for all SAT conditions
idxFastAccu = true(size(spkCorr,1),1);
idxFast = ismember(spkCorr.satCondition,'Fast');
idxAccu = ismember(spkCorr.satCondition,'Accurate');

idxErrNonErr = true(size(spkCorr,1),1);
idxError = spkCorr.isSefErrorUnit == 1;
idxNonError = spkCorr.isSefErrorUnit == 0;

sefUnitTypes = {'Error & non-Error', 'Error', 'Non-Error'};
conds = {'FAST & ACCURATE','FAST','ACCURATE'};
outcomes = {'ALL', 'Correct', 'ErrorChoice', 'ErrorTiming'};
rhoCondTbl = table();
rhoValTbl = table();
% rho itself
minMaxRho = minmax(spkCorr.rhoRaw');
rhoBins = -0.8:0.05:0.8; rhoBinEdges = [-0.75 rhoBins+0.05];

for ut = 1:numel(sefUnitTypes)
    sefUnitType = sefUnitTypes{ut};    
    if strcmp(sefUnitType,'Error')
        utIdx = idxError;
        sefUnitTypeStr = 'Error Units';
        oFilePrefix = 'Error_units_';
        yMaxs = [140 70 50 50]; % after seeing the plot
    elseif strcmp(sefUnitType,'Non-Error')
        utIdx = idxNonError;
        sefUnitTypeStr = 'Non-Error Units';
        oFilePrefix = 'Non_Error_units_';
        yMaxs = [140 70 50 50]; % after seeing the plot
    else
        utIdx = idxErrNonErr;
        sefUnitTypeStr = 'All Units';
        oFilePrefix = 'All_units_';
        yMaxs = [280 130 80 80]; % after seeing the plot
    end


for c = 1:numel(conds)
    cond = conds{c};    
    if strcmp(cond,'FAST')
        condIdx = idxFast;
        condStr = 'FAST';
        oFile = 'rhoDistrib_FAST.pdf';
        yMaxs = [140 70 50 50]; % after seeing the plot
    elseif strcmp(cond,'ACCURATE')
        condIdx = idxAccu;
        condStr = 'ACCURATE';
        oFile = 'rhoDistrib_ACCURATE.pdf';
        yMaxs = [140 70 50 50]; % after seeing the plot
    else
        condIdx = idxFastAccu;
        condStr = 'FAST & ACCURATE';
        oFile = 'rhoDistrib_FAST_ACCURATE.pdf';
        yMaxs = [280 130 80 80]; % after seeing the plot
    end
    idxPos = spkCorr.rhoRaw >=0 & condIdx & utIdx;
    idxNeg = spkCorr.rhoRaw < 0 & condIdx & utIdx;
    idxPosEst40 = spkCorr.rhoEstRaw_nTrials_40 >= 0 & condIdx & utIdx;
    idxNegEst40 = spkCorr.rhoEstRaw_nTrials_40 < 0 & condIdx & utIdx;
  
    % plot for each cond
    H_plots = getPlotHandles(); % 8
    nPlot = 0;
    useRhoCols = {'rhoRaw','rhoEstRaw_nTrials_40'}; 
    for ur = 1:numel(useRhoCols)
        rhoCol = useRhoCols{ur};
        rhoOrEstRho = 'Rho distribution for ';
        obsOrEstim = 'Observed';
        if contains(rhoCol,'Est')
           obsOrEstim = 'Estimated';         
        end
        
        idxPos = spkCorr.(rhoCol) >=0 & condIdx;
        idxNeg = spkCorr.(rhoCol) < 0 & condIdx;
        
        for oc = 1:numel(outcomes)
            outcome = outcomes{oc};
            idx = ismember(spkCorr.outcome,outcome);
            if sum(idx) == 0
                idx = true(size(spkCorr,1),1);
            end
            titleStr = upper(outcome);
            
            nPlot = nPlot + 1;
            axes(H_plots(nPlot));
            [ot] = plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNeg & idx),(idxPos & idx),yMaxs(oc));
            h_t = title(titleStr);
            if oc == 1
                % add observed or Estimated
                xlim = get(gca,'XLim');
                pos = get(h_t,'Position');
                text(xlim(1),pos(2),obsOrEstim,'FontWeight','bold','FontSize',14,'VerticalAlignment','bottom');
            end
            tempCondTbl = table();
            tempCondTbl.condition = {cond};
            tempCondTbl.outcome = {outcome};
            tempCondTbl.obsOrEstim = {obsOrEstim};
            rhoValTbl = [rhoValTbl;[tempCondTbl ot]];
            rhoCondTbl = [rhoCondTbl; tempCondTbl];
        end
    end
    %side annotation
    h=axes(gcf,'Position',[0.02,0.5,0.0001,0.0001]);
    txtStr = sprintf('[%s]  %s  -  %s',sefUnitTypeStr, rhoOrEstRho ,condStr);
    h_txt = text(0,0,txtStr,'Rotation',90,'FontWeight','bold','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle');
    saveFigPdf([oFilePrefix oFile]);
end

end


%%

for c = 1:numel(conds)
    cond = conds{c};
    
    if strcmp(cond,'FAST')
        condIdx = idxFast;
        oFile = 'rhoDistrib_FAST.pdf';
        yMaxs = [140 70 50 50]; % after seeing the plot
    elseif strcmp(cond,'ACCURATE')
        condIdx = idxAccu;
        oFile = 'rhoDistrib_ACCURATE.pdf';
        yMaxs = [140 70 50 50]; % after seeing the plot
    else
        condIdx = idxFastAccu;
        oFile = 'rhoDistrib_FAST_ACCURATE.pdf';
        yMaxs = [280 130 80 80]; % after seeing the plot
    end
    
    idxPos = spkCorr.rhoRaw >=0 & condIdx;
    idxNeg = spkCorr.rhoRaw < 0 & condIdx;
    idxPosEst40 = spkCorr.rhoEstRaw_nTrials_40 >= 0 & condIdx;
    idxNegEst40 = spkCorr.rhoEstRaw_nTrials_40 < 0 & condIdx;
    
    
    
    
    H_plots = getPlotHandles(); % 8
    nPlot = 0;
    % for all
    rho = spkCorr.rhoRaw;
    absRho = abs(spkCorr.rhoRaw);
    sqrRho = spkCorr.rhoRaw.^2;
    
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,idxNeg,idxPos,yMaxs(1))
    title({'Dist. Rho values';'ALL'})
    % for correct
    idx = ismember(spkCorr.outcome,'Correct');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNeg & idx),(idxPos & idx),yMaxs(2))
    title('CORRECT')
    % for Error choice
    idx = ismember(spkCorr.outcome,'ErrorChoice');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNeg & idx),(idxPos & idx),yMaxs(3))
    title('ERROR CHOICE')
   % for Error timing
    idx = ismember(spkCorr.outcome,'ErrorTiming');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNeg & idx),(idxPos & idx),yMaxs(4))
    title('ERROR TIMING')
    
    % Estimated Rho for all
    rho = spkCorr.rhoEstRaw_nTrials_40;
    absRho = abs(spkCorr.rhoEstRaw_nTrials_40);
    sqrRho = spkCorr.rhoEstRaw_nTrials_40.^2;

    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,idxNegEst40,idxPosEst40,yMaxs(1))
    title({'Dist. Estimated Rho [random 40] values';'ALL'})
    
    % Estimated Rho for correct
    idx = ismember(spkCorr.outcome,'Correct');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNegEst40 & idx),(idxPosEst40 & idx),yMaxs(2))
    title('CORRECT')
    % Estimated Rho for Error choice
    idx = ismember(spkCorr.outcome,'ErrorChoice');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNegEst40 & idx),(idxPosEst40 & idx),yMaxs(3))
    title('ERROR CHOICE')
    % Estimated Rho for Error timing
    idx = ismember(spkCorr.outcome,'ErrorTiming');
    nPlot = nPlot + 1;
    axes(H_plots(nPlot));
    plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,(idxNegEst40 & idx),(idxPosEst40 & idx),yMaxs(4))
    title('ERROR TIMING')
    
    h=axes(gcf,'Position',[0.02,0.5,0.0001,0.0001]);
    h_txt = text(0,0,cond,'Rotation',90,'FontWeight','bold','FontSize',16);
    
    saveFigPdf(oFile);
end

%% nTrials distribution
minMaxTrials = minmax(spkCorr.nTrials');
nTrialBins = 15:30:645; nTrialBinEdges = [0 nTrialBins+15];
[rhoNegNTrialC]=histcounts(spkCorr.nTrials(idxNeg),nTrialBinEdges);
[rhoPosNTrialC]=histcounts(spkCorr.nTrials(idxPos),nTrialBinEdges);

figure
bar(nTrialBins,[rhoNegNTrialC;rhoPosNTrialC]')
xlabel('nTrials')
ylabel('frequency')
legend({'neg. rho','pos. rho'},'Location','northeast')
title('Distribution of num. trials for Negative and positive rho values')

%%
function [outTbl] = plotRhoDist(rhoBins,rhoBinEdges,rho,absRho,sqrRho,idxNeg,idxPos,yLimMax)
[rhoNegC]=histcounts(rho(idxNeg),rhoBinEdges);
[rhoPosC]=histcounts(rho(idxPos),rhoBinEdges);
idxNegPos = idxNeg | idxPos;

bar(rhoBins,[rhoNegC;rhoPosC]')
xlabel('rho')
ylabel('frequency')
legend({'neg. rho','pos. rho'},'Location','northeast')

outTbl = table();
outTbl.countAbsolute = sum(idxNegPos);
outTbl.absRho_mean	= nanmean(absRho(idxNegPos));
outTbl.absRho_mean_sqr = nanmean(sqrRho(idxNegPos));
outTbl.countNegtiv	= sum(idxNeg);
outTbl.negtivRho_mean = nanmean(rho(idxNeg)); 
outTbl.negtivRho_mean_sqr = nanmean(sqrRho(idxNeg));
outTbl.countPostiv	= sum(idxPos);
outTbl.postivRho_mean = nanmean(rho(idxPos));
outTbl.postivRho_mean_sqr = nanmean(sqrRho(idxPos));


anno = { sprintf('Num. values\n  Rho-Neg = %d\n  Rho-Pos = %d',outTbl.countNegtiv,outTbl.countPostiv)
    '-----------------------------------'
    'Mean Rho Values:'
    sprintf('Abs = %0.3f, sqrd = %0.3f', outTbl.absRho_mean,outTbl.absRho_mean_sqr)
    sprintf('Neg = %0.3f, sqrd = %0.3f', outTbl.negtivRho_mean,outTbl.negtivRho_mean_sqr)
    sprintf('Pos = %0.3f, sqrd = %0.3f', outTbl.postivRho_mean,outTbl.postivRho_mean_sqr)
    };
set(gca,'YLim',[0,yLimMax]); % after seeing the plots
yLim = get(gca,'YLim');
h_txt = text(rhoBins(1),yLim(2)*0.95,anno,'VerticalAlignment', 'top');
end

%%
function [H_plots] = getPlotHandles()

nRows = 2;
nCols = 4;
    H_Figure = newFigure();
    pltH = 0.4;
    pltW = 0.20;
    gutter = 0.04;
    offsetsX = 0.06:(pltW+gutter):1-gutter; % for 4 column-starts
    offsetsY = 0.95-pltH:-(pltH+gutter+0.01):gutter; % for 3 row-starts
    H_plots = [];
    pltCount = 0;
    for ro = 1:nRows
        pos(2) = offsetsY(ro);
        for col = 1:nCols
            pos(1) = offsetsX(col);
            pos(3:4) = [pltW pltH];
            
            pltCount = pltCount + 1;
            H_plots(pltCount) = axes('parent',H_Figure,'position',pos,'box','on', 'layer','top','Tag',sprintf('H_plot%d',pltCount)); %#ok<AGROW>
        end
    end

end

%%
function [colNames] = getColNamesToUse()
colNames = {
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