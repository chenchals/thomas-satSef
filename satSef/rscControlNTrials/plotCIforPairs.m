%%
spkCorr = load('newRscWithSubSampling.mat');
spkCorr = spkCorr.spkCorr;

rscData = spkCorr(:,{'condition','rhoRaw_150ms',...
                     'nTrials','nTrials4SubSample',...
                     'normalCI95_150ms','percentileCI95_150ms',...
                     'rhoRawWithinNormalCI','rhoRawWithinPercentileCI'}...
                 );
 colNames = rscData.Properties.VariableNames;
 rscData.Properties.VariableNames = regexprep(colNames,'_150ms','');

conditions = {'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
              'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
    };
rscData.absRho = abs(rscData.rhoRaw);
rscData.satCondition = regexprep(rscData.condition,{'Correct','Error.*'},{'',''});
rscData.outcome = regexprep(rscData.condition,{'Fast','Accurate'},{'',''});

rscData = sortrows(rscData,'rhoRaw');

%%
plotRscWithCI(rscData,'normal','rscForPairsWithCI')
ppretty([10,16])

%% 
rscDataStats = grpstats(rscData(:,{'satCondition','outcome','absRho','nTrials'}),{'satCondition','outcome'},{'mean','sem'});
rscDataStats = sortrows(rscDataStats,{'outcome','satCondition'});
rscDataStats.Properties.RowNames = {};
% display 3 groups of 2 bars each
% data = [1 2;3 4;5 6]
%rsc groups
% rscDataStats =
% 
%   6×7 table
% 
%     satCondition       outcome       GroupCount    mean_absRho    sem_absRho    mean_nTrials    sem_nTrials
%     ____________    _____________    __________    ___________    __________    ____________    ___________
% 
%      'Accurate'     'Correct'           258         0.082981      0.0042224        469.67         6.2587   
%      'Fast'         'Correct'           258         0.082178      0.0046837        409.88         5.9226   
%      'Accurate'     'ErrorChoice'       258           0.1174      0.0054808        79.558         2.1948   
%      'Fast'         'ErrorChoice'       258          0.12106      0.0068834        169.77         3.1139   
%      'Accurate'     'ErrorTiming'       258          0.10216      0.0054097        153.83         2.8459   
%      'Fast'         'ErrorTiming'       258           0.1621      0.0082477        37.279         1.0289   
% quick:
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
sat = {'Accurate','Fast'};
meanRho = [0.0830    0.0822;    0.1174    0.1211;    0.1022    0.1621];
semRho = [0.0042    0.0047;    0.0055    0.0069;    0.0054    0.0082];
meanNTrials = [469.6667  409.8760;   79.5581  169.7713;  153.8295   37.2791];
semNtrials =[6.2587    5.9226;    2.1948    3.1139;    2.8459    1.0289];
figure
groupedBarPlot(meanRho,semRho,'Spike count corr (Rsc)')
ylabel('Spike Count Corr (Rsc)')
figure
groupedBarPlot(meanNTrials,semNtrials,'Trial count (Rsc)')
ylabel('Num. Trials')

%% trial distribution for pairs for Fast and Accurate Error timing:
idx = ismember(rscData.condition,'FastErrorTiming');
FastErrorTimTrialCouns = rscData.nTrials(idx);
idx = ismember(rscData.condition,'AccurateErrorTiming');
AccurateErrorTimTrialCouns = rscData.nTrials(idx);
Z=[FastErrorTimTrialCouns;AccurateErrorTimTrialCouns];
accClr = [1 0.2 0.2];
fasClr = [0.2 1.0 0.2];


bc = 0:10:250;
histogram(FastErrorTimTrialCouns,bc,'FaceColor',fasClr)
hold on
histogram(AccurateErrorTimTrialCouns,bc,'FaceColor',accClr)
ppretty([4,5])





%%
function [] = groupedBarPlot(meanVals, semVals, titleStr)
accClr = [1 0.2 0.2];
fasClr = [0.2 1.0 0.2];

hb = bar(meanVals,'grouped');
set(hb(1),'FaceColor',accClr);
set(hb(2),'FaceColor',fasClr);

hold on
ngroups = 3;
nbars = 2;
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x,meanVals(:,i), semVals(:,i), 'k', 'linestyle', 'none');
end
title(titleStr)
ppretty([2.5,2.5])

end

%% Plot Rsc with confidence intervals for each pair by condition
function [] = plotRscWithCI(rscData, normalOrPercentileCI, figName)
    conditions = {'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
                  'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
        };

    f = newFigure();
    for cn = 1:numel(conditions)
        condition = conditions{cn};
        idx = ismember(rscData.condition,condition);
        condTbl = rscData(idx,:);

        rho = condTbl.rhoRaw;
        %% lines
        x = 1:size(condTbl,1);
        if strcmpi(normalOrPercentileCI,'normal')
            ciStr = 'normal';
            ci = cell2mat(condTbl.normalCI95)';
            inOutCI = condTbl.rhoRawWithinNormalCI;
        elseif strcmpi(normalOrPercentileCI,'percentile')

            ciStr = 'percentile';
            ci = cell2mat(condTbl.percentileCI95)';
            inOutCI = condTbl.rhoRawWithinPercentileCI;
        end
        figFile = [figName '_' ciStr '.pdf'];

        withinCI = sum(inOutCI == 1);
        outsideCI = sum(inOutCI == 0);

        maxCi = max(ci(:));
        xx = [x;x;nan(1,numel(x))];
        yy = [ci;nan(1,numel(x))];    
        xx=xx(:);
        yy = yy(:);  

        subplot(6,1,cn)
        hl = line(xx,yy,'Color',[1 0.2 0.2],'Linewidth',1.5);
        hold on
        scatter(x,rho,30,'.','b')
        str = sprintf('%s within CI = %d, outside CI = %d',condition,withinCI,outsideCI);
        text(10,maxCi*0.9,str)
        ylabel('Rsc')
        xlabel('Pair count (sorted by Rsc)')
        xlim([0 260])
        legend({'Rsc CI 95%','Rsc'},'Location','southeast')

        if (cn == 1)
            title(sprintf('Spike count correlation with Confidence Interval (%s)',ciStr))
        end
    end
saveFigPdf(figFile);
end
