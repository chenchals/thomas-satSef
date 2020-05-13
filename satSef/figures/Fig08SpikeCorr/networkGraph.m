%%draw a dummy graph
objGr = dummyGraph()

%% Read and trim data for network graph
spkCorr = load('dataProcessed/satSefPaper/rscSubSampl1K_PostSaccade.mat');
% significant spk corr values
spkCorr = spkCorr.spkCorr;
% recode error neurons
spkCorr.isSefErrorUnit = abs(spkCorr.X_errGrade) > 1 | abs(spkCorr.X_rewGrade) > 1;
% recode sat and outcome
spkCorr.satCondition = regexprep(spkCorr.condition,{'Correct','Error.*'},{'',''});
spkCorr.outcome = regexprep(spkCorr.condition,{'Fast','Accurate'},{'',''});
% spike count corr sign
spkCorr.plusRho(spkCorr.rhoRaw > 0) = 1;
spkCorr.minusRho(spkCorr.rhoRaw < 0) = 1;

%% Filter
% 
% Considering ONLY SEF error (choice + timing) neurons, and ONLY those with
% stat signif correlation values, for each session draw a diagram like that
% below.  Significant correlations will be referred to as ?connections?  
idxPlus = spkCorr.isSefErrorUnit == 1 ...
      & ismember(spkCorr.outcome,{'ErrorChoice','ErrorTiming'}) ...
      & spkCorr.signifRaw_05 == 1; 
spkCorr = spkCorr(idxPlus,{'X_monkey','X_sess','satCondition','outcome','X_unitNum','Y_area','plusRho','minusRho'});
spkCorr = sortrows(spkCorr,{'X_unitNum','satCondition','outcome','Y_area','plusRho','minusRho'});

%% For each SEF error unit, find how many FEF or SC units they are connected to
[emptyGraph] = getEmptyGraph();
satConds = {'Fast','Accurate'};
outcomes = {'ErrorChoice_ErrorTiming','ErrorChoice','ErrorTiming'};
outStruct = struct();
for sc = 1:numel(satConds)
    satCondition = satConds{sc};
    for oc = 1:numel(outcomes)
        connGraph = emptyGraph;
        outcome = outcomes{oc};
        idxPlus = ismember(spkCorr.satCondition,satCondition) ...
            & contains(outcome,spkCorr.outcome);
        countTbl = spkCorr(idxPlus,{'X_unitNum','Y_area','plusRho','minusRho'});
        sumTbl = grpstats(countTbl,{'X_unitNum','Y_area'},{'sum'});
        
        sefErrUnits = unique(sumTbl.X_unitNum,'stable');
        for un = 1:numel(sefErrUnits)
            unitNum = sefErrUnits(un);
            unitTbl = sumTbl(sumTbl.X_unitNum == unitNum,:);
            if size(unitTbl,1) == 2 % unit is connected to both FEF and SC
                srcNode = 'Err2FEF_SC';
                % FEF target
                idxPlus = connGraph.findedge(srcNode,'FEF_Plus');
                idxMinus = connGraph.findedge(srcNode,'FEF_Minus');
                idxSum = find(ismember(unitTbl.Y_area,'FEF'));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);
                % SC target
                idxPlus = connGraph.findedge(srcNode,'SC_Plus');
                idxMinus = connGraph.findedge(srcNode,'SC_Minus');
                idxSum = find(ismember(unitTbl.Y_area,'SC'));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);                  
            else % unit is connected to either FEF or SC
                targNodePre = unitTbl.Y_area{1};
                srcNode = ['Err2' targNodePre];
                % FEF or SC target
                idxPlus = connGraph.findedge(srcNode,[targNodePre '_Plus']);
                idxMinus = connGraph.findedge(srcNode,[targNodePre '_Minus']);
                idxSum = find(ismember(unitTbl.Y_area,targNodePre));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);
             end
        end
        outStruct.(satCondition).(outcome).connGraph = connGraph;       
        outStruct.(satCondition).(outcome).countTbl = countTbl;
    end
end

% 






%%
function [emptyGraph] = getEmptyGraph()
nodeNames = {
     'Err2FEF'   
     'Err2FEF_SC'
     'Err2SC'    
     'FEF_Plus'  
     'FEF_Minus'
     'SC_Plus' 
     'SC_Minus' 
    };
nodePositions = [
    1,3
    1,2
    1,1
    2,2.6
    2,2.4
    2,1.6
    2,1.4];
nodesTbl = cell2table(nodeNames,'VariableNames',{'nodeName'});
nodesTbl.nodeId = (1:size(nodesTbl,1))';
srcNodeIds = [1 1 2 2 2 2 3 3];% lookup into node Ids
targNodeIds = [4 5 4 5 6 7 6 7];% lookup into nodeIds
weights = [0 0 0 0 0 0 0 0]; % counts 
tempArr = {
    'Err2FEF','FEF_Plus','error unit to FEF: positive Rsc';
    'Err2FEF','FEF_Minus','error unit to FEF: negative Rsc';
    'Err2FEF_SC','FEF_Plus','error unit to FEF and SC: FEF-positive Rsc';
    'Err2FEF_SC','FEF_Minus','error unit to FEF and SC: FEF-negative Rsc';
    'Err2FEF_SC','SC_Plus','error unit to FEF and SC: SC-positive Rsc';
    'Err2FEF_SC','SC_Minus','error unit to FEF and SC: SC-negative Rsc';
    'Err2SC','SC_Plus','error unit to SC: positive Rsc';
    'Err2SC','SC_Minus','error unit to SC: negative Rsc';
    };
emptyGraph = graph(srcNodeIds,targNodeIds,weights,nodesTbl.nodeName);
% linestyle
emptyGraph.Edges.lineStyle = repmat({'-'},numel(weights),1);
idxMinus = contains(emptyGraph.Edges.EndNodes(:,2),'Minus');
emptyGraph.Edges.lineStyle(idxMinus) = repmat({'--'},sum(idxMinus),1);
% edgecolor
emptyGraph.Edges.edgeColor = repmat({[0 0 1]},numel(weights),1);
emptyGraph.Edges.edgeColor(idxMinus) = repmat({[1 0 1]},sum(idxMinus),1);
% comments
emptyGraph.Edges.Comment = tempArr(:,3);
% node positions on graph
emptyGraph.Nodes.xPos = nodePositions(:,1);
emptyGraph.Nodes.yPos = nodePositions(:,2);

% 
%              EndNodes              Weight                      Comment                   
%     ___________________________    ______    ____________________________________________
% 
%     'Err2FEF'       'FEF_Plus'       0       'error unit to FEF: positive Rsc'           
%     'Err2FEF'       'FEF_Minus'      0       'error unit to FEF: negative Rsc'           
%     'Err2FEF_SC'    'FEF_Plus'       0       'error unit to FEF and SC: FEF-positive Rsc'
%     'Err2FEF_SC'    'FEF_Minus'      0       'error unit to FEF and SC: FEF-negative Rsc'
%     'Err2FEF_SC'    'SC_Plus'        0       'error unit to FEF and SC: SC-positive Rsc' 
%     'Err2FEF_SC'    'SC_Minus'       0       'error unit to FEF and SC: SC-negative Rsc' 
%     'Err2SC'        'SC_Plus'        0       'error unit to SC: positive Rsc'            
%     'Err2SC'        'SC_Minus'       0       'error unit to SC: negative Rsc'            
% 
%%
end

%%
function [objGr] = dummyGraph()
objGr = getEmptyGraph();
objGr.Edges.Weight = randi(12,[size(objGr.Edges,1),1]);
%lw = scaleVector(objGr.Edges.Weight,1,10);
lw = objGr.Edges.Weight;
plot(objGr,'XData',objGr.Nodes.xPos,'YData',objGr.Nodes.yPos,...
    'LineStyle',objGr.Edges.lineStyle,...
    'LineWidth',lw,...
    'EdgeColor',cell2mat(objGr.Edges.edgeColor),...
    'EdgeLabel',objGr.Edges.Weight,...
    'Interpreter','None');
end

function [vecScaled] = scaleVector(vec,minLim, maxLim)   
     minVec = min(vec);
     maxVec = max(vec);
     vecScaled = (vec - minVec)./(maxVec-minVec);
     vecScaled = (maxLim-minLim) * vecScaled + minLim;
end