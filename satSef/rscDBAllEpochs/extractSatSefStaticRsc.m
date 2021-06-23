% Extract static Rsc data for all SEF/FEF-SC pairs
% see also CREATESATSEFRSCDB

oFilename = 'dataProcessed/satSefPaper/analysis/spkCorr/summary/SAT_SEF_StaticRscAllPairs.mat';
areaPairs = {
    'SEF-FEF'    
    'SEF-SC'     
    'FEF-SC'     
    'SEF-SEF' 
    'FEF-FEF'    
    'SC-SC'      
    };
corrMatDirs = strcat('dataProcessed/satSefPaper/analysis/spkCorr/spkCorr_',areaPairs,'/mat');
corrDatFields = {
    % cellPairInfo fields to extract
    'Pair_UID'        
    'X_monkey'        
    'X_sessNum'       
    'X_sess'          
    'X_unitNum'       
    'Y_unitNum'       
    'X_unit'          
    'Y_unit'          
    'X_area'          
    'Y_area'          
    'X_visGrade'      
    'Y_visGrade'      
    'X_visField'      
    'Y_visField'      
    'X_visType'       
    'Y_visType'       
    'X_moveGrade'     
    'Y_moveGrade'     
    'X_moveField'     
    'Y_moveField'     
    'X_errGrade'      
    'Y_errGrade' 
    'X_isErrGrade'      
    'Y_isErrGrade' 
    'X_errField'      
    'Y_errField'      
    'X_rewGrade'      
    'Y_rewGrade'      
    'X_isRewGrade'      
    'Y_isRewGrade'      
    'X_taskType'      
    'X_Hemi'          
    'Y_Hemi'          
    'X_Grid'          
    'Y_Grid'          
    'X_GridAP_ML'     
    'Y_GridAP_ML'     
    'X_Depth'         
    'Y_Depth'         
    'X_Depth0'        
    'Y_Depth0'        
    'X_newDepth'      
    'Y_newDepth'      
    'XY_Dist'         
    'isOnSameChannel' 
    % spkCorr fields to extract
    'condition'          
    'alignedName'        
    'alignedEvent'       
    'alignedTimeWin'     
    'trialNosByCondition'
    'xBaselineMeanStd' 
    'yBaselineMeanStd'
    'rho_pval_win_150ms' 
    'critRho10'          
    'critRho05'          
    'critRho01' 
    'rho_pval_static_150ms'    
    'critRho10_Z_baseline'
    'critRho05_Z_baseline'
    'critRho01_Z_baseline'
    'rho_pval_static_Z_baseline_150ms'
    };
t = regexp(corrDatFields,'rho_pval_static_(\d*)ms$','tokens');
staticWinSizes = sort(cellfun(@(x) str2double(x{1}),t(~cellfun(@isempty,t))));
spkCorrStatic = struct();
nThreads = 0;
p = gcp('nocreate');
if numel(p) == 1
    nThreads = p.NumWorkers;
end
tic
for d = 1:numel(corrMatDirs)
    areaPair = areaPairs{d};
    fprintf('Doing pairs for %s...',areaPair);
    areaPairField = strrep(areaPair,'-','_');
    srcFiles = dir([corrMatDirs{d},'/spkCorr_PAIR_*.mat']);
    srcFiles = strcat({srcFiles.folder}','/',{srcFiles.name}');
    nPairs = numel(srcFiles);
    outPairs = struct();
    parfor (p = 1:nPairs, nThreads)
    %for p = 1:nPairs
        pDat = table();
        srcFile = srcFiles{p};
        pDatStruct = load(srcFile);
        [~,fn,ext] = fileparts(srcFile);
        srcFile = [fn ext]; 
        nRows = size(pDatStruct.spkCorr,1);
        temp = [repmat(pDatStruct.cellPairInfo,nRows,1) pDatStruct.spkCorr];
        tempFns = temp.Properties.VariableNames';
        pDat.srcFile = repmat({srcFile},nRows,1);
        pDat.pairAreas = repmat({areaPair},nRows,1);
        pDat.nTrials = cellfun(@(x) numel(x), temp.trialNosByCondition);
        % temp.trialNosByCondition = [];
        fieldIdx = cell2mat(cellfun(@(x) find(strcmp(tempFns,x)),corrDatFields,'UniformOutput',false));
        pDat = [pDat,temp(:,corrDatFields)];
        pDat.XY_Dist = cell2mat(pDat.XY_Dist);
        % split rho, pval
        for sW = 1:numel(staticWinSizes)
            staticWin = staticWinSizes(sW);
            swSuffix = num2str(staticWin,'_%dms');
            % Use rho/pval from raw counts            
            swRhoPval = pDat.(['rho_pval_static' swSuffix]);
            [rho,pval,sig05,sig01] = cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),swRhoPval);            
            [pDat.(['rhoRaw' swSuffix]),...
            pDat.(['pvalRaw' swSuffix]),...
            pDat.(['signifRaw_05' swSuffix]),...
            pDat.(['signifRaw_01' swSuffix])] = deal(rho,pval,sig05,sig01);
            % Use rhp/pval from Z-scored (with mean and std from) Baseline period
            swRhoPval = pDat.(['rho_pval_static_Z_baseline' swSuffix]);
            [rho,pval,sig05,sig01] = cellfun(@(x) deal(x(1),x(2),x(2)<=0.05,x(2)<=0.01),swRhoPval);            
            [pDat.(['rhoZBaseline' swSuffix]),...
            pDat.(['pvalZBaseline' swSuffix]),...
            pDat.(['signifZBaseline_05' swSuffix]),...
            pDat.(['signifZBaseline_01' swSuffix])] = deal(rho,pval,sig05,sig01);
        end
                
        outPairs(p).pDat = pDat;
    end
    
    spkCorrStatic.(areaPairField) = vertcat(outPairs.pDat); 
    fprintf('Done %.3f sec.\n',toc)
end
size(spkCorrStatic);
[d,~,~]=fileparts(oFilename);
if ~exist(d,'dir')
    mkdir(d);
end
save(oFilename,'-v7.3','-struct','spkCorrStatic');
toc

