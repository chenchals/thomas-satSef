function [anovaResults] = satAnova(valsGroupsTbl,anoveModelName,doMultCompareFlag,alpha)
%SATANOVA Do a multiway anova for the given table of inputs
%   Column1 = Y-Values numeric
%   Columns (2 to end) = groups/factors over which the anonan is run
% ***Terminology: group == factor
% ***Note: dont use "factor" *** --> Its is a matlab funtion for primes  doc FACTOR
% see also ANOVAN, MULTCOMPARE

%anoveModelName = 'linear'; %'interaction'
%doMultCompareFlag = true;

anovaDisplay = 'off'; %'on'

yVals = valsGroupsTbl{:,1};
assert(isnumeric(yVals),'Y-Values must be numeric');

groupNames = valsGroupsTbl.Properties.VariableNames;
groupNames = groupNames(2:end);

groups = cell(1,numel(groupNames));
for ii = 1:numel(groupNames)
    groups{ii} = valsGroupsTbl.(groupNames{ii});
end

% anovaTbl header:
% {'Source' 'Sum Sq.' 'd.f.' 'Singular?' 'Mean Sq.' 'F'  'Prob>F'}
anovaResults = struct();
anovaTblVarNames = {'Source', 'SumSq' 'df' 'IsSingular' 'MeanSq' 'F'  'ProbGtF'};

[~,temp,anovaStats] = anovan(yVals,groups,'model',anoveModelName,'varnames',groupNames, 'display',anovaDisplay);
anovaTbl = cell2table(temp(2:end,:),'VariableNames',anovaTblVarNames);
% add '*' p(F >= .05) and '**' p(F >=  .01)
idx = find(~ismember(anovaTbl.Source,{'Error','Total'}));
for jj = 1:numel(idx)
    str = '';
    probGtF = anovaTbl.ProbGtF{jj};
    if probGtF <= 0.01
        str = '**';
    elseif probGtF<=0.05
        str = '*';
    end
    anovaTbl.signifStr{jj} = str;
end
anovaResults.anovaTbl = anovaTbl;
% Compare results for different LEVELS *WITHIN* each group/Factor independently 
% for bonferroni use 'CType', ... see doc multcompare
if (doMultCompareFlag)
    for gr = 1:numel(groupNames)
        [temp,~,~,grpNames] = multcompare(anovaStats,'Dimension',gr,'Alpha',alpha);
        anovaResults.(groupNames{gr}) = annotateMultcompareResults(temp,grpNames);
    end
    
    % Compare results for different LEVEL combinations *ACROSS* each group/Factor independently
    nWays = 2;
    n2GrpComparisions = combnk(1:numel(groupNames),nWays);
    
    for jj = 1:size(n2GrpComparisions,1)
        idx = n2GrpComparisions(jj,:);
        fn = char(join(groupNames(idx),'_'));
        [temp,~,~,grpNames] = multcompare(anovaStats,'Dimension',idx,'Alpha',alpha);
        anovaResults.(fn) = annotateMultcompareResults(temp,grpNames);
    end
end

%% If there are 3 groups also do a 3 way comparision
if numel(groupNames) > 2
    nWays = 3;
    n3GrpComparisions = combnk(1:numel(groupNames),nWays);

    for jj = 1:size(n3GrpComparisions,1)
        idx = n3GrpComparisions(jj,:);
        fn = char(join(groupNames(idx),'_'));
        [temp,~,~,grpNames] = multcompare(anovaStats,'Dimension',idx,'Alpha',alpha);
        anovaResults.(fn) = annotateMultcompareResults(temp,grpNames);
    end
end

end % fxn : satAnova()

function [resultsAsTbl] = annotateMultcompareResults(results,grpNames)
% see also multcompare

levelNamesNew = split(grpNames,{',','='});
levelNamesNew = levelNamesNew(:,2:2:end);
levelNamesNew = arrayfun(@(x) char(join(levelNamesNew(x,:),'_')),1:size(levelNamesNew,1),'UniformOutput',false)';

resultsAsTbl = table();
resultsAsTbl.levelName1 = arrayfun(@(x) levelNamesNew(x),results(:,1),'UniformOutput',false);
resultsAsTbl.levelName2 = arrayfun(@(x) levelNamesNew(x),results(:,2),'UniformOutput',false);
resultsAsTbl = [resultsAsTbl array2table(results(:,3:end),'VariableNames',{'loCI95','meanDiff','hiCI95','pval_H0'})];
resultsAsTbl.isSignif05 = results(:,end)<=0.05;
resultsAsTbl.isSignif01 = results(:,end)<=0.01;
resultsAsTbl.signifStr=repmat({''},size(resultsAsTbl,1),1);
resultsAsTbl.signifStr(resultsAsTbl.isSignif05) = repmat({'*'},sum(resultsAsTbl.isSignif05),1);
resultsAsTbl.signifStr(resultsAsTbl.isSignif01) = repmat({'**'},sum(resultsAsTbl.isSignif01),1);
resultsAsTbl.group1 = results(:,1);
resultsAsTbl.group2 = results(:,2);

end % util : annotateMultcompareResults()
