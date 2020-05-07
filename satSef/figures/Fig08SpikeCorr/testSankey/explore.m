
% see: https://www.r-graph-gallery.com/323-sankey-diagram-with-the-networkd3-library.html
% data from the url 'https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json'
% txt is a struct
txt = webread('https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json');

% Read Rsc Data for PostSaccade from excel file
rscData = readtable('fig08_data.xlsx');


outcome = 'ErrorTiming';
idx = ismember(rscData.outcome,outcome);
myData = rscData(idx,:);
% code Rsc sign: 1 = negative; 2 = positive
myData.codedSign(myData.rscObserved < 0) = 1;
myData.codedSign(myData.rscObserved > 0) = 2;
% separate FAST and ACCURATE
fastTbl = myData(ismember(myData.satCondition,'Fast'),{'PairUid','isSefErrorUnit','codedSign','rscObserved'});
accuTbl = myData(ismember(myData.satCondition,'Accurate'),{'PairUid','isSefErrorUnit','codedSign','rscObserved'});

temp = outerjoin(accuTbl,fastTbl,'LeftKeys',{'PairUid'},'RightKeys',{'PairUid'});
fastAccuTbl = table();
fastAccuTbl.PairUid = temp.PairUid_accuTbl;
% code SEF node name
idx = temp.isSefErrorUnit_accuTbl==1;
fastAccuTbl.sefErrorTye = repmat({''},size(temp,1),1);
fastAccuTbl.sefErrorTye(idx) = repmat({'SEF-Error'},sum(idx),1);
fastAccuTbl.sefErrorTye(~idx) = repmat({'SEF-Other'},sum(~idx),1);
% code Pair's Other node for  Sat cond and sign
fastAccuTbl.AccuRhoSign = repmat({''},size(temp,1),1);
idx = temp.codedSign_accuTbl == 1; % negative Rsc
fastAccuTbl.AccuRho = temp.rscObserved_accuTbl;
fastAccuTbl.AccuRhoSign(idx) = repmat({'AccuMinus'},sum(idx),1);
fastAccuTbl.AccuRhoSign(~idx) = repmat({'AccuPlus'},sum(~idx),1);
% do for fast
fastAccuTbl.FastRhoSign = repmat({''},size(temp,1),1);
idx = temp.codedSign_fastTbl == 1; % negative Rsc
fastAccuTbl.FastRho = temp.rscObserved_fastTbl;
fastAccuTbl.FastRhoSign(idx) = repmat({'FastMinus'},sum(idx),1);
fastAccuTbl.FastRhoSign(~idx) = repmat({'FastPlus'},sum(~idx),1);

fastAccuTbl.codedSignFast = temp.codedSign_fastTbl;
fastAccuTbl.codedSignAccu = temp.codedSign_accuTbl;
fastAccuTbl.codedSignSum2Sort = fastAccuTbl.codedSignAccu + fastAccuTbl.codedSignFast;
% sort the table
fastAccuTbl = sortrows(fastAccuTbl,{'sefErrorTye','codedSignSum2Sort','codedSignAccu','codedSignFast','PairUid'});

%%
% nodes: (Left) Accu Positive, Accu Negative, 
%       (Middle) SEF ErrorNeurons, SEF OtherNeurons, 
%       (Right) Fast Positive, Fast Negative
sank = struct();
nodeNames = unique([fastAccuTbl.AccuRhoSign;fastAccuTbl.sefErrorTye;fastAccuTbl.FastRhoSign],'stable')';
sank.nodes = cell2struct(nodeNames,'name',1);
sourceTarget = {{'AccuMinus','SEF-Error'}
                {'AccuMinus','SEF-Other'}
                {'AccuPlus','SEF-Error'}
                {'AccuPlus','SEF-Other'}
                {'SEF-Error','FastMinus'}
                {'SEF-Error','FastPlus'}
                {'SEF-Other','FastMinus'}
                {'SEF-Other','FastPlus'}
                };
links = struct();            
for l = 1:numel(sourceTarget)
    src = sourceTarget{l}{1};
    targ = sourceTarget{l}{2};
    source = find(strcmp(nodeNames,src)) - 1;
    target = find(strcmp(nodeNames,targ)) - 1;
    if contains(src,'Accu')
        idx = ismember(fastAccuTbl.AccuRhoSign,src) & ismember(fastAccuTbl.sefErrorTye,targ);
        rhoVal = mean(fastAccuTbl.AccuRho(idx));
        counts =  sum(idx);
    elseif contains(targ,'Fast')
        idx = ismember(fastAccuTbl.FastRhoSign,targ) & ismember(fastAccuTbl.sefErrorTye,src);
        rhoVal = mean(fastAccuTbl.FastRho(idx));
        counts =  sum(idx);
    end
    links(l,1).source = source;
    links(l,1).target = target;
    links(l,1).counts = counts; 
    links(l,1).rhoValMean = rhoVal; 
    
end
sank.links=links;

% write json file
fid = fopen('satSef/figures/Fig08SpikeCorr/testSankey/connErrorOther.json','w');
fwrite(fid,jsonencode(sank));
fclose(fid);


%%
