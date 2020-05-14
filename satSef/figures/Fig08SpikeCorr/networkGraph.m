%%draw a dummy graph
%  objGr = getEmptyGraph();
%  objGr.Edges.Weight = randi(20,[size(objGr.Edges,1),1]);
% % 
%  h = plotGraph(objGr);

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
idxErr = spkCorr.isSefErrorUnit == 1 ...
      & ismember(spkCorr.outcome,{'ErrorChoice','ErrorTiming'}) ...
      & spkCorr.signifRaw_05 == 1; 
spkCorr = spkCorr(idxErr,{'X_monkey','X_sess','satCondition','outcome','X_unitNum','Y_area','plusRho','minusRho'});
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
        idxSat = ismember(spkCorr.satCondition,satCondition);
        if strcmp(outcome,'ErrorChoice_ErrorTiming')
            idxOutcome = ones(size(spkCorr,1),1);
        else
            idxOutcome = ismember(spkCorr.outcome,outcome);
        end
        countTbl = spkCorr(idxSat & idxOutcome,{'X_unitNum','Y_area','plusRho','minusRho'});
        sumTbl = grpstats(countTbl,{'X_unitNum','Y_area'},{'sum'});
        
        sefErrUnits = unique(sumTbl.X_unitNum,'stable');
        for un = 1:numel(sefErrUnits)
            unitNum = sefErrUnits(un);
            unitTbl = sumTbl(sumTbl.X_unitNum == unitNum,:);
            if size(unitTbl,1) == 2 % unit is connected to both FEF and SC
                srcNodeP = 'SEF2FEF_SC_P';
                srcNodeM = 'SEF2FEF_SC_M';
                % FEF target
                idxPlus = connGraph.findedge(srcNodeP,'FEF_Plus');
                idxMinus = connGraph.findedge(srcNodeM,'FEF_Minus');
                idxSum = find(ismember(unitTbl.Y_area,'FEF'));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);
                % SC target
                idxPlus = connGraph.findedge(srcNodeP,'SC_Plus');
                idxMinus = connGraph.findedge(srcNodeM,'SC_Minus');
                idxSum = find(ismember(unitTbl.Y_area,'SC'));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);                  
            else % unit is connected to either FEF or SC
                targNodePre = unitTbl.Y_area{1};
                srcNodeP = ['SEF2' targNodePre '_P'];
                srcNodeM = ['SEF2' targNodePre '_M'];
                % FEF or SC target
                idxPlus = connGraph.findedge(srcNodeP,[targNodePre '_Plus']);
                idxMinus = connGraph.findedge(srcNodeM,[targNodePre '_Minus']);
                idxSum = find(ismember(unitTbl.Y_area,targNodePre));
                connGraph.Edges.Weight(idxPlus) = connGraph.Edges.Weight(idxPlus) + unitTbl.sum_plusRho(idxSum);
                connGraph.Edges.Weight(idxMinus) = connGraph.Edges.Weight(idxMinus) + unitTbl.sum_minusRho(idxSum);
             end
        end
        outStruct.(satCondition).(outcome).connGraph = connGraph;       
        outStruct.(satCondition).(outcome).countTbl = countTbl;
    end
end
%% plot network graph for [all, ErrorChoice, ErrorTiming] for fast and accurate
H_axes = figTemplate(2,3);
plotNo = 0;
for sc = 1:numel(satConds)
    satCondition = satConds{sc};
    for oc = 1:numel(outcomes)
        outcome = outcomes{oc};
        objGraph = outStruct.(satCondition).(outcome).connGraph;
        plotNo = plotNo + 1;
        axes(H_axes(plotNo));
        plotGraph(objGraph);
        title([satCondition ' - ' outcome],'FontSize',12,'FontWeight','bold','Interpreter','none');
    end
end
fn ='Da_Eu_NetworkPlot.pdf';
print(fn,'-fillpage','-dpdf','-painters')


%%
function [hA] = figTemplate(ros,cols)
    figure;
    set(gcf,'Color',[1 1 1],'Position',[100 100 1400 800],'PaperOrientation','landscape');
    hA = tight_subplot(ros, cols, [.05 .05],[0.05 0.05],[0.05 0.05]);
    for ii = 1:numel(hA)
        axes(hA(ii))
        title(sprintf('Plot %d',ii));
    end
end

%%
function [emptyGraph] = getEmptyGraph()
nodeNames = {
     'SEF2FEF_P'   
     'SEF2FEF_M'      
     'SEF2FEF_SC_P'
     'SEF2FEF_SC_M'
     'SEF2SC_P'    
     'SEF2SC_M'    
     'FEF_Plus'  
     'FEF_Minus'
     'SC_Plus' 
     'SC_Minus' 
    };
nodePositions = [
    1,3.1
    1,3.0
    1,2.1
    1,2.0
    1,1.1
    1,1.0
    2,2.6
    2,2.5
    2,1.6
    2,1.5];
nodesTbl = cell2table(nodeNames,'VariableNames',{'nodeName'});
nodesTbl.nodeId = (1:size(nodesTbl,1))';
srcNodeIds = [1 2 3 4 3 4 5 6];% lookup into node Ids
targNodeIds = [7 8 7 8 9 10 9 10];% lookup into nodeIds
weights = [0 0 0 0 0 0 0 0]; % counts 
tempArr = {
    'SEF2FEF_P','FEF_Plus','SEF error unit to FEF: positive Rsc';
    'SEF2FEF_M','FEF_Minus','SEF error unit to FEF: negative Rsc';
    'SEF2FEF_SC_P','FEF_Plus','SEF error unit to FEF and SC: FEF-positive Rsc';
    'SEF2FEF_SC_M','FEF_Minus','SEF error unit to FEF and SC: FEF-negative Rsc';
    'SEF2FEF_SC_P','SC_Plus','SEF error unit to FEF and SC: SC-positive Rsc';
    'SEF2FEF_SC_M','SC_Minus','SEF error unit to FEF and SC: SC-negative Rsc';
    'SEF2SC_P','SC_Plus','SEF error unit to SC: positive Rsc';
    'SEF2SC_M','SC_Minus','SEF error unit to SC: negative Rsc';
    };
emptyGraph = graph(srcNodeIds,targNodeIds,weights,nodesTbl.nodeName);
% linestyle
emptyGraph.Edges.lineStyle = repmat({'-'},numel(weights),1);
idxMinus = contains(emptyGraph.Edges.EndNodes(:,2),'Minus');
emptyGraph.Edges.lineStyle(idxMinus) = repmat({':'},sum(idxMinus),1);
% edgecolor
emptyGraph.Edges.edgeColor = repmat({[0 0 1]},numel(weights),1);
emptyGraph.Edges.edgeColor(idxMinus) = repmat({[1 0 1]},sum(idxMinus),1);
% comments
emptyGraph.Edges.Comment = tempArr(:,3);
% node positions on graph
emptyGraph.Nodes.xPos = nodePositions(:,1);
emptyGraph.Nodes.yPos = nodePositions(:,2);

% % emptyGraph = 
% % 
% %   graph with properties:
% % 
% %     Edges: [8×5 table]
% %     Nodes: [10×3 table]
% % 
% % emptyGraph.Nodes
% % 
% % ans =
% % 
% %   10×3 table
% % 
% %          Name         xPos    yPos
% %     ______________    ____    ____
% % 
% %     'SEF2FEF_P'        1      3.1 
% %     'SEF2FEF_M'        1        3 
% %     'SEF2FEF_SC_P'     1      2.1 
% %     'SEF2FEF_SC_M'     1        2 
% %     'SEF2SC_P'         1      1.1 
% %     'SEF2SC_M'         1        1 
% %     'FEF_Plus'         2      2.6 
% %     'FEF_Minus'        2      2.5 
% %     'SC_Plus'          2      1.6 
% %     'SC_Minus'         2      1.5 
% % 
% % emptyGraph.Edges
% % 
% % ans =
% % 
% %   8×5 table
% % 
% %               EndNodes               Weight    lineStyle     edgeColor                          Comment                     
% %     _____________________________    ______    _________    ____________    ________________________________________________
% % 
% %     'SEF2FEF_P'       'FEF_Plus'       0          '-'       [1×3 double]    'SEF error unit to FEF: positive Rsc'           
% %     'SEF2FEF_M'       'FEF_Minus'      0          ':'       [1×3 double]    'SEF error unit to FEF: negative Rsc'           
% %     'SEF2FEF_SC_P'    'FEF_Plus'       0          '-'       [1×3 double]    'SEF error unit to FEF and SC: FEF-positive Rsc'
% %     'SEF2FEF_SC_P'    'SC_Plus'        0          '-'       [1×3 double]    'SEF error unit to FEF and SC: FEF-negative Rsc'
% %     'SEF2FEF_SC_M'    'FEF_Minus'      0          ':'       [1×3 double]    'SEF error unit to FEF and SC: SC-positive Rsc' 
% %     'SEF2FEF_SC_M'    'SC_Minus'       0          ':'       [1×3 double]    'SEF error unit to FEF and SC: SC-negative Rsc' 
% %     'SEF2SC_P'        'SC_Plus'        0          '-'       [1×3 double]    'SEF error unit to SC: positive Rsc'            
% %     'SEF2SC_M'        'SC_Minus'       0          ':'       [1×3 double]    'SEF error unit to SC: negative Rsc'   
% %
end

%%
function [h_graph] = plotGraph(objGr)
%lw = scaleVector(objGr.Edges.Weight,1,10);
% remoce all edges with 0 weight (ie no connections)
inValidEdgeIdx = find(objGr.Edges.Weight == 0);
objGr = objGr.rmedge(inValidEdgeIdx); %#ok<FNDSB>
lw = objGr.Edges.Weight;
h_graph = plot(objGr,'XData',objGr.Nodes.xPos,'YData',objGr.Nodes.yPos,...
    'LineStyle',objGr.Edges.lineStyle,...
    'LineWidth',lw,...
    'EdgeColor',cell2mat(objGr.Edges.edgeColor),...
    'EdgeLabel',objGr.Edges.Weight,...
    'Interpreter','None');
set(h_graph,'NodeLabel',{});
%set(h_graph,'Marker','s','MarkerSize',20,'NodeColor',[0.5 0.5 0.5])
bgColor = get(gcf,'Color');
set(get(h_graph,'Parent'),'XColor',bgColor,'YColor',bgColor)
set(get(h_graph,'Parent'),'XTick',[],'YTick',[]);
% annotation of nodes
text(0.85,3.05,'SEF to FEF','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','top');
text(0.85,2.05,'SEF to FEF & SC','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','top');
text(0.85,1.05,'SEF to SC','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','top');
text(2.15,2.55,'FEF','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(2.15,1.55,'SC','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','bottom');

end

function [vecScaled] = scaleVector(vec,minLim, maxLim)   
     minVec = min(vec);
     maxVec = max(vec);
     vecScaled = (vec - minVec)./(maxVec-minVec);
     vecScaled = (maxLim-minLim) * vecScaled + minLim;
end