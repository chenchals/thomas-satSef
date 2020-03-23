function [satSefPairCellInfoDB,satSefPairSummary] = createSatSefCellPairsInfoDB()
% 
% Needs 'dataNeurophys_SAT.mat' created by Thomas(10/29/2019)
%       and '[monk]_SAT_colorRecode.xlsx' that are color recoded versions
%       of Rich's SAT summary excel files
% Criteria for Pairs:
%      1. There has to be atleast 2 units in the session
%      2. The data for all units and their criteria is in file
%      'dataNeurophys_SAT.mat'  
% 
% After running this script, 
%       run createTrialTypesEventTimesDB 
% see also: CREATESATSEFTRIALTYPESEVENTTIMESDB, PARSESATEXCEL
% Modifications:
% 10/03/2019 : Removed the criteria for filtering for only vis/mov units
% 10/29/2019 : Using new file 'dataNeurophys_SAT.mat' for ninfo and nstats
%              instead of 'ninfo_nstats_SAT.mat'
% 03/23/2020 : New script from
%       scripts/pairAnalysis/createDataset/createJpsthPairCellInfoDB.m
%              Converted to a funtion to conserve memory for temp vars
%

%%
% baseSatSefDir All sat-sef-paer anaysis (new)
baseSatSefDir = 'dataProcessed/satSefPaper';
xlsxDir = [baseSatSefDir,'/excel'];
datasetDir = [baseSatSefDir,'/dataset'];
% if true then: goodVMIdx = abs(visGrade) > 1 | abs(moveGrade) > 1
% else use all units from cInfo 
% cInfo is from unitInfoStats_file
useVizMovUnits = false;
unitInfoStats_file = fullfile(datasetDir, 'dataNeurophys_SAT.mat');
% useNSEFNUnits : For sat sef Paper remove NSEFN units. 
useNSEFNUnits = false;

% translated and original datafiles for sessions
matAndPlxFiles = fullfile(datasetDir,'SessionFiles_SAT.mat');
% output files for SatSef Paper
satSefPairCellInfoDBFile = fullfile(datasetDir,'SAT_SEF_PAIR_CellInfoDB.mat');
satSefPairSummaryFile = fullfile(datasetDir,'SAT_SEF_PAIR_Summary.mat');
satSefPairSummaryFileCsv = fullfile(datasetDir,'SAT_SEF_PAIR_Summary.csv');

% Use only Darwin and Euler
monkNameMap = containers.Map({'D','E'},{'Darwin','Euler'});

%% Add matlab datafile and plx datafile  to the info table
sessMatPlxFiles = load(matAndPlxFiles);
sessMatPlxFiles = sessMatPlxFiles.sessionFiles;

%% Add information of grid location and hemifield from recoded excel files
excelInfos = cellfun(@(x) parseSatExcel(fullfile(xlsxDir,[x '_SAT_colorRecode.xlsx'])), monkNameMap.values, 'UniformOutput', false);
excelInfos = vertcat(excelInfos{:});

%% Process info for all monks 
db = load(unitInfoStats_file); % contains ninfo and nstats structs
cInfo = db.unitInfo;

% convert grid location from ap-ml to 
% +(plus) for anterior -(minus) for posterior
% +(plus) for medial -(minus) for lateral
fx_GridConvert = @(gloc) fliplr(eval([ '[' regexprep(fliplr(gloc),{'a','p','m','l'},{'+','-','+','-'}) ']']));
gridLoc = regexprep(lower(excelInfos.Grid),{'-','nan',' '},{''});
excelInfos.GridAP_ML = cellfun(@(x) fx_GridConvert(x),gridLoc,'UniformOutput',false);

newDepth = str2double(excelInfos.Depth);
d0Idx = isnan(str2double(excelInfos.Depth)) & ~isnan(str2double(excelInfos.Depth0));
newDepth(d0Idx) = str2double(excelInfos.Depth0(d0Idx));
excelInfos.newDepth = newDepth;

% inner join with information from excel recoded
cInfo = innerjoin(cInfo,excelInfos,'LeftKeys',{'sess','unit'},...
                  'RightKeys',{'MatSessionName','Unit'},...
                  'RightVariables',{'MatSessionName','Hemi','Grid','GridAP_ML',...
                  'Depth','Depth0','newDepth','SessionNotes'});
              
%% JPSTH summary --> how many cell pairs for each session?
JpsthPairSummary=table();
JpsthPairSummary.sess = unique(cInfo.sess);
JpsthPairSummary.nUnits = cell2mat(cellfun(@(x) sum(contains(cInfo.sess,x)),JpsthPairSummary.sess,'UniformOutput',false));
temp = innerjoin(sessMatPlxFiles,JpsthPairSummary);
JpsthPairSummary.matDatafile = temp.matDatafile;
JpsthPairSummary.plxDatafile = temp.plxDatafile;

if useVizMovUnits
    goodVMIdx = abs(cInfo.visGrade) > 1 | abs(cInfo.moveGrade) > 1; %#ok<UNRCH>
    goodUnits = cInfo(goodVMIdx,:);
else
    goodUnits = cInfo;
end
% filter out NSEFN units from goodUnits
if ~useNSEFNUnits
    goodUnits(strcmp(goodUnits.area,'NSEFN'),:) = [];
end

cellsBySession = arrayfun(@(x) find(contains(goodUnits.sess,x)), JpsthPairSummary.sess, 'UniformOutput',false);
varsForPairs = cInfo.Properties.VariableNames;
nextPairId = 0;
JpsthPairCellInfoDB = table();

%% Ensure the following is true for area1 vs area2
% for cross area pairs: always have first SEF on x-axis, then FEF on
% x-axis, then SC on x-axis; NSEFN is always on y-axis
% Retain NSEFN (not-SEF_neuron) so, for future, the pair IDs correspond
if useNSEFNUnits
    pairXYarea = {
        {'SEF' 'SEF'}
        {'SEF' 'FEF'}
        {'SEF' 'SC'}
        {'SEF' 'NSEFN'}
        {'FEF' 'FEF'}
        {'FEF' 'SC'}
        {'FEF' 'NSEFN'}
        {'SC' 'SC'}
        {'SC' 'NSEFN'}
        {'NSEFN' 'NSEFN'}
        };
else
    pairXYarea = {
        {'SEF' 'SEF'}
        {'SEF' 'FEF'}
        {'SEF' 'SC'}
        {'FEF' 'FEF'}
        {'FEF' 'SC'}
        {'SC' 'SC'}
        };
end
% concated strings
pairXYarea = cellfun(@(x) [x{:}],pairXYarea,'UniformOutput',false);
          
for s=1:numel(cellsBySession)
    res = goodUnits(cellsBySession{s},:);
    session = res.sess{1}; 
    fprintf('Processing session [%s]\n',session);
    tIdx = contains(JpsthPairSummary.sess,session);
    if size(res,1) <= 1
        JpsthPairSummary.nNSEFNUnits(tIdx) = 0;
        JpsthPairSummary.nSEFUnits(tIdx) = size(res(strcmp(res.area,'SEF'),1),1);
        JpsthPairSummary.nFEFUnits(tIdx) = size(res(strcmp(res.area,'FEF'),1),1);
        JpsthPairSummary.nSCUnits(tIdx) = size(res(strcmp(res.area,'SC'),1),1);
        JpsthPairSummary.nCellsForJpsth(tIdx) = 0;
        JpsthPairSummary.nPairsJpsth(tIdx) = 0;
        JpsthPairSummary.nSEF_SEF(tIdx) = 0;
        JpsthPairSummary.nSEF_FEF(tIdx) = 0;
        JpsthPairSummary.nSEF_SC(tIdx) = 0;
        JpsthPairSummary.nFEF_FEF(tIdx) = 0;
        JpsthPairSummary.nFEF_SC(tIdx) = 0;
        JpsthPairSummary.nSC_SC(tIdx) = 0;       
        continue;
    elseif size(res,1) > 1 % we have more than 1 unit
        result.CellInfoTable = goodUnits(cellsBySession{s},:);
        sessName = JpsthPairSummary.sess{tIdx};
        matDatafile = JpsthPairSummary.matDatafile{tIdx};
        plxDatafile = JpsthPairSummary.plxDatafile{tIdx};        
        
        pairRowIds = sortrows(combnk(1: size(result.CellInfoTable,1), 2),[1 2]);
        nPairs = size(pairRowIds,1);
        pairs = table();
        pairs.Pair_UID = cellstr(num2str(((1:nPairs)+ nextPairId)','PAIR_%04d'));
        
        % NSEFN units per session
        JpsthPairSummary.nNSEFNUnits(tIdx) = JpsthPairSummary.nUnits(tIdx) - size(result.CellInfoTable,1);
        JpsthPairSummary.nSEFUnits(tIdx) = size(res(strcmp(res.area,'SEF'),1),1);
        JpsthPairSummary.nFEFUnits(tIdx) = size(res(strcmp(res.area,'FEF'),1),1);
        JpsthPairSummary.nSCUnits(tIdx) = size(res(strcmp(res.area,'SC'),1),1);
        JpsthPairSummary.nCellsForJpsth(tIdx) = size(result.CellInfoTable,1);
        JpsthPairSummary.nPairsJpsth(tIdx) = size(pairs,1);
               
        % add flag to check if we need to swap rowIds so as to get X-Unit
        % and Y-Uint too be congruent with pairXYAreas defined above
        % the pairRowIds(:,3) is the swar flag, swap col1 and 2 for that
        % row
        pairRowIds(:,3) = arrayfun(@(x) ...
            sum(strcmp(pairXYarea,[result.CellInfoTable.area{pairRowIds(x,:)}]))== 0,...
            1:nPairs)';
        pairRowIds(:,4:5) = pairRowIds(:,1:2);
        swapCols = find(pairRowIds(:,3));
        for zz = 1:numel(swapCols)
            pairRowIds(swapCols(zz),1:2) = fliplr(pairRowIds(swapCols(zz),1:2));
        end
        for v = 1:numel(varsForPairs)
            cName = varsForPairs{v};
            pairs.(['X_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,1));
            pairs.(['Y_' cName]) = result.CellInfoTable.(cName)(pairRowIds(:,2));
            if strcmp(cName,'errGrade')
               pairs.('X_isErrGrade') = abs(pairs.(['X_' cName]))>=2;
               pairs.('Y_isErrGrade') = abs(pairs.(['Y_' cName]))>=2;
            end
            if strcmp(cName,'rewGrade')
               pairs.('X_isRewGrade') = abs(pairs.(['X_' cName]))>=2;
               pairs.('Y_isRewGrade') = abs(pairs.(['Y_' cName]))>=2;
            end
            
        end
        XY_Dist = arrayfun(@(x)...
            getElectrodeDistance(pairs.X_GridAP_ML{x},pairs.Y_GridAP_ML{x},...
                                pairs.X_newDepth(x),pairs.Y_newDepth(x)),...
                                (1:size(pairs,1))','UniformOutput',false);
        % set cross area distances to be NaN as the distance computation is
        % invalid across chambers
        XY_Dist(cellfun(@(a1,a2) ~strcmp(a1,a2),pairs.X_area,pairs.Y_area)) = {NaN};
        pairs.XY_Dist = XY_Dist;
        
        pairs.isOnSameChannel = arrayfun(@(x) ...
            isequal(regexp(pairs.X_unit{x},'(\d{1,2})','tokens'),...
            regexp(pairs.Y_unit{x},'(\d{1,2})','tokens')), (1:size(pairs))');
        
        pairs.matDatafile = repmat({matDatafile},nPairs,1);          
        pairs.plxDatafile = repmat({plxDatafile},nPairs,1);          
        JpsthPairCellInfoDB = [JpsthPairCellInfoDB;pairs]; %#ok<AGROW>
        
        % add to summary table
        JpsthPairSummary.nSEF_SEF(tIdx) = fx_nArea1Area2Pairs(pairs,'SEF','SEF');
        JpsthPairSummary.nSEF_FEF(tIdx) = fx_nArea1Area2Pairs(pairs,'SEF','FEF');
        JpsthPairSummary.nSEF_SC(tIdx) = fx_nArea1Area2Pairs(pairs,'SEF','SC');
        JpsthPairSummary.nFEF_FEF(tIdx) = fx_nArea1Area2Pairs(pairs,'FEF','FEF');
        JpsthPairSummary.nFEF_SC(tIdx) = fx_nArea1Area2Pairs(pairs,'FEF','SC');
        JpsthPairSummary.nSC_SC(tIdx) = fx_nArea1Area2Pairs(pairs,'SC','SC');      

        
        nextPairId = nextPairId + nPairs;
    end
    result.PairInfoTable = pairs;
end
%% Save output files
satSefPairCellInfoDB = JpsthPairCellInfoDB;
satSefPairSummary= JpsthPairSummary;
save(satSefPairCellInfoDBFile,'satSefPairCellInfoDB');
save(satSefPairSummaryFile, 'satSefPairSummary');
writetable(satSefPairSummary,satSefPairSummaryFileCsv);

end

%% Sub functions
function [nPairs] = fx_nArea1Area2Pairs(pairs,xArea,yArea)
  nPairs = sum(ismember(pairs.X_area,xArea) & ismember(pairs.Y_area,yArea));
end

function [eDist] = getElectrodeDistance(xGridLoc,yGridLoc,xDepth,yDepth)
% grid locs must be numeric [AP, ML]
% [+1 = 1a, -1 = 1p] [+1 = 1m, -1 = 1L]
% depth in microns --> convert to mm

    if isempty(xGridLoc) || isempty(yGridLoc) || sum(isnan([xDepth yDepth]))
        eDist = NaN;
    elseif isequal(xGridLoc,yGridLoc)
        eDist = 0;
    else
        pointA = [xGridLoc(:); xDepth/1000];
        pointB = [yGridLoc(:); yDepth/1000];
        eDist =  sqrt(sum((pointA-pointB).^2));
    end

end
