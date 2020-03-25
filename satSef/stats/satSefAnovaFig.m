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
%%
% compute anova for the following epochs
%spkCorr.rho = spkCorr.rhoZ;

% for our plots: 1 way anova for each condition is used
statsAnovaF1Outcome = doAnovaOutcome(spkCorr,'rhoZ');
statsAnova.AccurateOutcomeTbl = statsAnovaF1Outcome.Accurate_PostSaccade.outcome;
statsAnova.AccurateAnovaTbl = statsAnovaF1Outcome.Accurate_PostSaccade.anovaTbl;
statsAnova.AccurateStatsTbl = statsAnovaF1Outcome.Accurate_PostSaccade.statsTbl;

statsAnova.FastOutcomeTbl = statsAnovaF1Outcome.Fast_PostSaccade.outcome;
statsAnova.FastAnovaTbl = statsAnovaF1Outcome.Fast_PostSaccade.anovaTbl;
statsAnova.FastStatsTbl = statsAnovaF1Outcome.Fast_PostSaccade.statsTbl;
clearvars statsAnovaF1Outcome useCols

%%

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
function [anovaConditionByOutcome] = doAnovaConditionByOutcome(spkCorr)
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
        anovaTbl.rho = spkCorr{idx,'rho'};
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