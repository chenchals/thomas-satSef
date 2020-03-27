%load r_sc data
% Use onkly SEF_FEF and SEF_SC pairs
load('dataProcessed/satSefPaper/analysis/spkCorr/summary/SAT_SEF_StaticRscAllPairs.mat', 'SEF_FEF','SEF_SC');
useCols = {'pairAreas','condition','alignedName','rhoRaw_150ms','rhoZBaseline_150ms'};
% useCols
spkCorr = [SEF_FEF(:,useCols); SEF_SC(:,useCols)];
% recode columns for anova
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
spkCorr.epoch = spkCorr.alignedName;
% take absolute value of correlation
spkCorr.rhoR = abs(spkCorr.rhoRaw_150ms);
spkCorr.rhoZ = abs(spkCorr.rhoZBaseline_150ms);

%% One-Way Anova foer outcomes for each SAT condition
[statsAnovaOutcome] = computeOneWay(spkCorr,'rhoR');

%% Two-way for SatCondition by Outcome
% main effect : Sat Condition
% Interaction term: condition by outcome: Specifically:
%   Accurate_Correct vs Fast_Correct, 
%   Accurate_ErroeChoice vs Fast_ErrorChoice,
%   Accurate_ErrorTiming vs Fast_ErrorTiming
statsAnovaF2CondByOutcome = doAnovaConditionByOutcome(spkCorr,'rhoZ');

%%

function [statsAnova] = computeOneWay(spkCorr, useRho)
%% compute 1-way for coutcome for each SAT condition separately 
statsAnovaF1Outcome = doAnovaOutcome(spkCorr,useRho);
% flatten tables - Accurate
statsAnova.AccurateOutcomeTbl = statsAnovaF1Outcome.Accurate_PostSaccade.outcome;
statsAnova.AccurateAnovaTbl = statsAnovaF1Outcome.Accurate_PostSaccade.anovaTbl;
statsAnova.AccurateStatsTbl = statsAnovaF1Outcome.Accurate_PostSaccade.statsTbl;
% flatten tables - Fast
statsAnova.FastOutcomeTbl = statsAnovaF1Outcome.Fast_PostSaccade.outcome;
statsAnova.FastAnovaTbl = statsAnovaF1Outcome.Fast_PostSaccade.anovaTbl;
statsAnova.FastStatsTbl = statsAnovaF1Outcome.Fast_PostSaccade.statsTbl;

% plot - outcomes for Fast Condition
meanRhoCol = ['mean_' useRho];
semRhoCol = ['sem_' useRho];
% find max value
f = statsAnova.FastStatsTbl;
a = statsAnova.AccurateStatsTbl;
[af] = [a.(meanRhoCol);f.(meanRhoCol)];
maxY = max(af)*1.2;
figure
fastClr = [0 0.7 0];
bar(categorical(f.outcome),f.(meanRhoCol),'FaceColor',fastClr,'EdgeColor','none')
hold on
ebF = errorbar(categorical(f.outcome),f.(meanRhoCol),f.(semRhoCol),'LineStyle','none','LineWidth',1,'Color','k');
ylabel('R_{sc} SEF-FEF/SC Pairs')
ylim([0 maxY])
ppretty(2)
title('Fast','Position',[0.5,maxY,0]);

% print - Fast - stats for replotting / significance
fastLines = [strsplit(evalc('disp(statsAnova.FastStatsTbl)'),'\n')'
         strsplit(evalc('disp(statsAnova.FastOutcomeTbl)'),'\n')'
         strsplit(evalc('disp(statsAnova.FastAnovaTbl)'),'\n')'];
fastLines = regexprep(fastLines,{'<strong>','</strong>'},{'',''});
fprintf('--------------FAST - STATS - LINES-------------------------\n')
fastLines
fprintf('-----------------------------------------------------------\n')

% plot - outcomes for Accurate condition
figure
accClr = [0.8 0 0];
p2 = bar(categorical(a.outcome),a.(meanRhoCol),'FaceColor',accClr,'EdgeColor','none')
hold on
ebA = errorbar(categorical(a.outcome),a.(meanRhoCol),a.(semRhoCol),'LineStyle','none','LineWidth',1,'Color','k');
ylabel('R_{sc} SEF-FEF/SC Pairs')
ylim([0 maxY])
ppretty(2)
title('Accurate','Position',[0.5,maxY,0]);

% print - Accurate - stats for replotting / significance
accuLines = [strsplit(evalc('disp(statsAnova.AccurateStatsTbl)'),'\n')'
         strsplit(evalc('disp(statsAnova.AccurateOutcomeTbl)'),'\n')'
         strsplit(evalc('disp(statsAnova.AccurateAnovaTbl)'),'\n')'];
accuLines = regexprep(accuLines,{'<strong>','</strong>'},{'',''});

fprintf('-----------ACCURATE - STATS - LINES------------------------\n')
accuLines
fprintf('-----------------------------------------------------------\n')

end




function [statsAnovaF1Outcome,statsAnovaF2ConditionByOutcome] = doAllAnova(spkCorr)
    statsAnovaF1Outcome = doAnovaOutcome(spkCorr);
    statsAnovaF2ConditionByOutcome = doAnovaConditionByOutcome(spkCorr);
end

function [anovaOutcome] = doAnovaOutcome(spkCorr,useRho)
%  Do 1-Factor anova for the following -
% .   Factor1: outcome
% The data colums for input anovaTbl from spkCorr data above:
%   yVals  : rho or rhoZ
%   outcome : Factor for 1-way anova values:
%             Correct|ErrorChoice|ErrorTiming

    satConditions = {'Accurate', 'Fast'};
    epochs = {'PostSaccade'};
    anovaOutcome = struct();
    for sc = 1:numel(satConditions)
        satCondition = satConditions{sc};
        idx = ismember(spkCorr.satCondition,satCondition);
        for ep = 1:numel(epochs)
            epoch = epochs{ep};
            idx = ismember(spkCorr.epoch,epoch) & idx;
            anovaTbl = table();
            anovaTbl.(useRho) = spkCorr{idx,useRho};
            anovaTbl.outcome = spkCorr{idx,'outcome'};
            % do anova
            anovaModel = 'interaction'; %  [linear|interaction]
            multiCompareFlag = 1; % [0|1]
            alpha = 0.05;
            fn = [satCondition '_' epoch];
            anovaOutcome.(fn) = satAnova(anovaTbl,anovaModel,multiCompareFlag,alpha);
            anovaOutcome.(fn).statsTbl = grpstats(anovaTbl,{'outcome'},{'mean','std','sem'});
        end
    end
end
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