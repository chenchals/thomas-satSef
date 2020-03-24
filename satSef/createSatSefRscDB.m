% Create spike count correlation data set for pairs of units recorded from
% same session for same as well as different areas
% see also CREATECORRSPKDATASET

areaPairs ={
    'SEF','FEF'
    'SEF','SC'
    'FEF','SC'
    'SEF','SEF'
    'FEF','FEF'
    'SC','SC'
    };
tic
for ii = 1:size(areaPairs,1)
    createRscDB(areaPairs{ii,1},areaPairs{ii,2});
    toc
end
toc

% Sub functions **********

function [] = createRscDB(area1,area2)
% Modified from thomas-jpsth/scripts/pairAnalysis/rscCreate/createSpkCorrDataset.m
% Create complete dataset for all the pairs matching the criteria for
% area1, area2
% Expects the following files in specified location:
% Inofrmation about all possible cell pairs:
%        'dataProcessed/satSefPaper/dataset/SAT_SEF_PAIR_CellInfoDB.mat'
% Trial types of all sessions (Accurate, Fast, Correct,...):
%        'dataProcessed/satSefPaper/dataset/SAT_SEF_TrialTypesDB.mat'
% Event times for all trials and all sessions:
%         'dataProcessed/satSefPaper/dataset/SAT_SEF_TrialEventTimesDB.mat'
% Get resptime and set it as SaccadePrimaryTempo
%         'dataProcessed/dataset/binfo_moves_SAT.mat'
% Spike time data for all units of all sessions:
%         'dataProcessed/dataset/spikes_SAT.mat'
    warning('off');
    monkIdsToDo = {'D','E'};
    %Options for Spk Corr computation
    rootAnalysisDir = 'dataProcessed/satSefPaper/analysis/spkCorr';
    datasetDir = 'dataProcessed/satSefPaper/dataset';
    % wavDir = 'dataProcessed/dataset/wavesNew';
    % Files for getting data to run JPSTH
    pairsFile = fullfile(datasetDir,'SAT_SEF_PAIR_CellInfoDB.mat');
    trialTypesFile = fullfile(datasetDir,'SAT_SEF_TrialTypesDB.mat');
    trialEventTimesFile = fullfile(datasetDir,'SAT_SEF_TrialEventTimesDB.mat');
    spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');

    % alignment:
    % Setup time windows for different event time alignment, the field names
    % SHALL correspond to column names for trialEventTimes below.
    alignNames = {'Baseline','Visual','PostSaccade','PostReward'};
    alignEvents = {'CueOn','CueOn','SaccadePrimary','RewardTime'};
    alignTimeWin = {[-600 100],[-200 400],[-100 500],[-200 700]};
    % first sort by accurate/fast?
    firstSortEvent = {'SaccadePrimary','SaccadePrimary','SaccadePrimary','SaccadePrimary'};
    %conditions
    conditionsTbl = table();
    conditionsTbl.conditions = {
        'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
        'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
        };

    resultsDir = fullfile(rootAnalysisDir,['spkCorr_' area1 '-' area2],'mat');
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end

    % Load data variable: JpsthPairsCellInfo
    cellPairs = load(pairsFile);
    cellPairs = cellPairs.satSefPairCellInfoDB;
    cellPairs = cellPairs(ismember([cellPairs.X_monkey],monkIdsToDo),:);
    % Load data variable: spike times
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    % Load data variable: TrialTypes
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    % Load data variable: TrialEventTimes
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;

    % Filter cell pairs for the areas of interest
    cellPairs = cellPairs(...
        ((strcmp(cellPairs.X_area,area1) & strcmp(cellPairs.Y_area,area2)) ...
        | (strcmp(cellPairs.X_area,area2) & strcmp(cellPairs.Y_area,area1))),...
        :);
    sessions = unique(cellPairs.X_sess);
    assert(isequal(sessions,unique(cellPairs.Y_sess)),'********Fatal: Error X-Unit sessions and Y-Unit sessions do not match');

    % For each JPSTH cell pair do SpikeCorr
    % see doc pctRunOnAll
    pctRunOnAll warning off;
    nPairs = size(cellPairs,1);
    nThreads = 0;
    p = gcp('nocreate');
    if numel(p) == 1
        nThreads = p.NumWorkers;
    end
    parfor (p = 1:nPairs,nThreads)
        cellPair = cellPairs(p,:); %#ok<*PFBNS>
        sess = cellPair.X_sess{1};
        % must be cell array of ntrials by 1
        xSpkTimes = spikesSat{cellPair.X_unitNum}';
        ySpkTimes = spikesSat{cellPair.Y_unitNum}';
        evntTimes = sessionEventTimes(strcmp(sessionEventTimes.session,sess),:);
        trialTypes = sessionTrialTypes(strcmp(sessionTrialTypes.session,sess),:);
        % No waveforms during analysis for SAT-SEF paper 
        spkCorr = getStaticRscForPair(cellPair,xSpkTimes,ySpkTimes,...
            evntTimes,trialTypes,conditionsTbl.conditions,...
            alignNames,alignEvents,alignTimeWin,firstSortEvent...
            );
        oFn = fullfile(resultsDir,['spkCorr_' cellPair.Pair_UID{1} '.mat']);
        saveSpkCorrData(oFn,spkCorr);
    end
end

function [] = saveSpkCorrData(oFn,varToSave)
    tic
    fprintf('Saving file : %s ...',oFn)
    tempConditions = varToSave;
    save(oFn,'-v7.3','-struct','tempConditions');
    fprintf('%d\n',toc);
end


function [outSpkCorr] = getStaticRscForPair(cellPair,xSpkTimes,ySpkTimes,evntTimes,trialTypes,...
    conditions,alignNames,alignEvents,alignTimeWin,firstSortEvent)
% Modified from thomas-jpsth/scripts/pairAnalysis/rscCreate/getSpkCorrForPair.m
% cellPair : the pair row in the table JPSTH_PAIRS_CellInfoDB.mat
% xSpkTimes : X-axis cell: spike times as cell array of nTrials
% ySpkTimes : Y-axis cell: spike times as cell array of nTrials
% evntTimes : event times for session from TrialEventTimesDB.mat
% trialTypes : trial types for session from TrialTypesDB.mat
% conditions : Trial type names: AccurateCorrect, FastCorrect... from
%              TrialTypesDB
%   Example:  conditions = {
%            'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
%            'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
%             };
% alignNames : Name of the aligned event window : Baseline, Visual,
%              postSaccade
% alignEvents : Name of the event to align on ... from
%               TrialEventTimesDB.mat
% alignTimeWin : A cell array of time windows for epochs corres[ponding to
%                align events

    % Static spk corr windows for computing spike corr
    % count spikes in this window : delta = 150 ms
    staticWinsAlignTs(1).Baseline = [-150 0];
    staticWinsAlignTs(1).Visual = [0 150];
    staticWinsAlignTs(1).PostSaccade = [0 150];
    staticWinsAlignTs(1).PostReward = [0 150];
 
    warning('off');
    % ignore processing if the sel. trials are below thisNum.
    % We need atleast 2 points to find Pearson's Correlation coefficient
    nTrialsThreshold = 2;

    units = struct();
    spkCorr = struct();

    XCellId = ['DSP' cellPair.X_unit{1}];
    YCellId = ['DSP' cellPair.Y_unit{1}];

    units.(XCellId) = xSpkTimes;
    units.(YCellId) = ySpkTimes;
    
    units.('X_trialMeanStd') = SpikeUtils.getTrialMeanStd(xSpkTimes);
    units.('Y_trialMeanStd') = SpikeUtils.getTrialMeanStd(ySpkTimes);
    
    baselineWinForMeanStd = [-600 -100];
    temp = SpikeUtils.alignSpikeTimes(xSpkTimes,evntTimes.CueOn{1},baselineWinForMeanStd);
    temp = SpikeUtils.rasters(temp,baselineWinForMeanStd);
    temp = temp.rasters;
    units.('X_baselineMeanStd') = arrayfun(@(x) [mean(temp(x,:)) std(temp(x,:))],(1:size(temp,1))','UniformOutput',false );
    temp = SpikeUtils.alignSpikeTimes(ySpkTimes,evntTimes.CueOn{1},baselineWinForMeanStd);
    temp = SpikeUtils.rasters(temp,baselineWinForMeanStd);
    temp = temp.rasters;
    units.('Y_baselineMeanStd') = arrayfun(@(x) [mean(temp(x,:)) std(temp(x,:))],(1:size(temp,1))','UniformOutput',false );
    clear temp;
    for cond = 1:numel(conditions)
        try
            % incase something breaks continue...
            condition = conditions{cond};
            selTrials = trialTypes.(condition){:};
            if isempty(selTrials)
                spkCorr.(condition) = [];
                continue;
            end

            % Mutually Exclusive trials for Choice/Timing Errors
            % If condition is *ChoiceErr or *TimingError ensure mutually
            % exclusive
            otherCondition = [];
            if contains(condition,'ChoiceError')
                otherCondition = regexprep(codition,'ChoiceError','TimingError');
            elseif contains(condition, 'TimingError')
                otherCondition = regexprep(codition,'TimingError','ChoiceError');
            end
            if ~isempty(otherCondition)
                selTrials(trialTypes.(otherCondition){:}) = 0;
            end

            % first check for nTrials
            if isempty(selTrials) || numel(selTrials) < nTrialsThreshold
                spkCorr.(condition) = [];
                continue;
            end

            % check if trials need to be dropped due to poor_isolation..
            trRem = cellPair.X_trRemSAT{1};
            if ~isempty(trRem)
                selTrials(trRem(1):trRem(2)) = 0;
            end
            trRem = cellPair.Y_trRemSAT{1};
            if ~isempty(trRem)
                selTrials(trRem(1):trRem(2)) = 0;
            end
            if isempty(selTrials)
                spkCorr.(condition) = [];
                continue;
            end
            % second check for nTrials
            selTrials = find(selTrials);
            if numel(selTrials) <= nTrialsThreshold
                spkCorr.(condition) = [];
                continue;
            end

            % Sort selected trials based on the event
            selTrlsTbl = table();
            selTrlsTbl.selTrials = selTrials;
            % for each aligned event
            opts = struct();
            for evId = 1:numel(alignEvents)
                alignedEvent = alignEvents{evId};
                alignedTimeWin = alignTimeWin{evId};
                alignedName = alignNames{evId};
                if isempty(firstSortEvent{evId})
                    firstSortName = [];
                    selTrlsTbl.firstSortTime = -inf(size(selTrlsTbl,1),1);
                else
                    firstSortName = firstSortEvent{evId};
                    firstSortTime = double(evntTimes.(firstSortName){1});
                    temp = firstSortTime(selTrlsTbl.selTrials);
                    temp(temp==0 | isnan(temp)) = -inf;
                    selTrlsTbl.firstSortTime = temp;
                end
                selTrlsTblSorted = sortrows(selTrlsTbl,{'firstSortTime'});
                selTrialsSorted = selTrlsTblSorted.selTrials;
                alignTime = evntTimes.CueOn{1};
                if ~strcmp(alignedEvent,'CueOn')
                    alignTime = alignTime + double(evntTimes.(alignedEvent){1}(:));
                end
                % Align Spike times and get rasters
                alignTime = alignTime(selTrialsSorted);
                XAligned = SpikeUtils.alignSpikeTimes(units.(XCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                YAligned = SpikeUtils.alignSpikeTimes(units.(YCellId)(selTrialsSorted),alignTime, alignedTimeWin);
                tempRast = SpikeUtils.rasters(XAligned,alignedTimeWin);
                XRasters = tempRast.rasters;
                tempRast = SpikeUtils.rasters(YAligned,alignedTimeWin);
                YRasters = tempRast.rasters;
                rasterBins = tempRast.rasterBins;                               
                clear tempRast;
                % Z scores
                % get Z score for given sample using means and stds specified for each trial
                % Z-Score using baseline mean and std
                fx_getZscore = @(rast,meanStd) cell2mat(arrayfun(@(x) (double(rast(x,:))- meanStd{1}(1))./meanStd{1}(2), (1:size(rast,1))','UniformOutput',false));
                xRasters_Z_baseline = fx_getZscore(XRasters,units.X_baselineMeanStd);
                yRasters_Z_baseline = fx_getZscore(YRasters,units.Y_baselineMeanStd);
                xRasters_Z_trial = fx_getZscore(XRasters,units.X_trialMeanStd);
                yRasters_Z_trial = fx_getZscore(YRasters,units.Y_trialMeanStd);

                % gather vars fo corrSpk_PAIR_xxxx.mat
                opts(evId,1).condition = {condition};
                opts(evId,1).alignedName = {alignedName};
                opts(evId,1).alignedEvent = {alignedEvent};
                opts(evId,1).alignedTimeWin = {alignedTimeWin};
                opts(evId,1).alignTime = {alignTime};
                opts(evId,1).firstSortByName = {firstSortName};
                opts(evId,1).firstSortByTime = selTrlsTblSorted.firstSortTime;
                opts(evId,1).trialNosByCondition = selTrialsSorted;
                opts(evId,1).xCellSpikeTimes = {XAligned}; %#ok<*AGROW>
                opts(evId,1).yCellSpikeTimes = {YAligned};
                opts(evId,1).xBaselineMeanStd = {cell2mat(units.X_baselineMeanStd)};
                opts(evId,1).yBaselineMeanStd = {cell2mat(units.Y_baselineMeanStd)};
                opts(evId,1).xTrialMeanStd = {cell2mat(units.X_trialMeanStd)};
                opts(evId,1).yTrialMeanStd = {cell2mat(units.Y_trialMeanStd)};
                opts(evId,1).rasterBins = {rasterBins};
                opts(evId,1).xRasters = {XRasters};
                opts(evId,1).yRasters = {YRasters};
                opts(evId,1).xRasters_Z_baseline = {xRasters_Z_baseline};
                opts(evId,1).yRasters_Z_baseline = {yRasters_Z_baseline};
                opts(evId,1).xRasters_Z_trial = {xRasters_Z_trial};
                opts(evId,1).yRasters_Z_trial = {yRasters_Z_trial};
                                
                % Compute spike corrs for static windows for each aligned event
                for sWin = 1:numel(staticWinsAlignTs)
                    staticWin = staticWinsAlignTs(sWin).(alignedName);
                    fieldSuffix = num2str(range(staticWin),'_%dms'); 
                    
                    % Static windows spike corr for Raw counts
                    opts(evId,1).(['rho_pval_win' fieldSuffix]) = {staticWin};
                    xSpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{XRasters},{staticWin},'UniformOutput',false);
                    ySpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{YRasters},{staticWin},'UniformOutput',false);
                    %rho_pval = {getCorrData(xSpkCounts{1},ySpkCounts{1},'Pearson')};
                    % we need the critical value of rho for .10, .05, .01 significance
                    [rho_pval1,opts(evId,1).critRho10,opts(evId,1).critRho05,opts(evId,1).critRho01] =...
                           getCorrData(xSpkCounts{1},ySpkCounts{1},'Pearson');
                    opts(evId,1).(['xSpkCount_win' fieldSuffix]) = xSpkCounts;
                    opts(evId,1).(['ySpkCount_win' fieldSuffix]) = ySpkCounts;
                    opts(evId,1).(['rho_pval_static' fieldSuffix]) = {rho_pval1};
                    
                    % Static windows spike corr for - Z-scored (using Baseline mean/std) count
                    xSpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{xRasters_Z_baseline},{staticWin},'UniformOutput',false);
                    ySpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{yRasters_Z_baseline},{staticWin},'UniformOutput',false);                     
                    
                    % We also need the critical valur of rho for .10, .05,
                    % .01 significance
                    %rho_pval = {getCorrData(xSpkCounts{1},ySpkCounts{1},'Pearson')};
                    [rho_pval2,opts(evId,1).critRho10_Z_baseline,opts(evId,1).critRho05_Z_baseline,opts(evId,1).critRho01_Z_baseline] = ...
                        getCorrData(xSpkCounts{1},ySpkCounts{1},'Pearson');
                                        
                    opts(evId,1).(['xSpkCount_win_Z_baseline' fieldSuffix]) = xSpkCounts;
                    opts(evId,1).(['ySpkCount_win_Z_baseline' fieldSuffix]) = ySpkCounts;
                    opts(evId,1).(['rho_pval_static_Z_baseline' fieldSuffix]) = {rho_pval2};
                    
                    % Static windows spike corr for - Z-scored (using Baseline mean/std) count
                    xSpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{xRasters_Z_trial},{staticWin},'UniformOutput',false);
                    ySpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                        {rasterBins},{yRasters_Z_trial},{staticWin},'UniformOutput',false);
                    rho_pval = {getCorrData(xSpkCounts{1},ySpkCounts{1},'Pearson')};

                    opts(evId,1).(['xSpkCount_win_Z_trial' fieldSuffix]) = xSpkCounts;
                    opts(evId,1).(['ySpkCount_win_Z_trial' fieldSuffix]) = ySpkCounts;
                    opts(evId,1).(['rho_pval_static_Z_trial' fieldSuffix]) = rho_pval;
                end
                                
            end % for alignEvents
            %tempSpkCorr.Properties.RowNames = alignNames;
            spkCorr.(condition) = struct2table(opts,'AsArray',true);
        catch mE
            getReport(mE)
            continue
        end
    end % for conditions
    % Save for the current pair
    fns = fieldnames(spkCorr);
    tempTbl = table();
    for fn = 1:numel(fns)
        if ~isempty(spkCorr.(fns{fn}))
            tempTbl = [tempTbl;spkCorr.(fns{fn})];
        end
    end
    % convert all doubles to single to save space
    fns = tempTbl.Properties.VariableNames;
    for ii = 1:numel(fns)
        fn = fns{ii};
        if iscell(tempTbl.(fn)) && isa(tempTbl.(fn){1},'double')
            tempTbl.(fn) = cellfun(@(x) single(x),tempTbl.(fn),'UniformOutput',false);
        elseif isa(tempTbl.(fn)(1),'double')
            tempTbl.(fn) = single(tempTbl.(fn));
        end
    end

    outSpkCorr.spkCorr = tempTbl;
    outSpkCorr.cellPairInfo = cellPair;
    outSpkCorr.GithubRef = 'https://github.com/chenchals/thomas-satSef.git';

end

function [rho_pval,critRho10,critRho05,critRho01] = getCorrData(xMat,yMat,corrMethodStr)
    % Get rho, pval from matlab corr function
    if strcmpi(corrMethodStr,'Pearson')
        corrMethod = 'Pearson';
    elseif strcmpi(corrMethodStr,'Spearman')
        corrMethod = 'Spearman';
    elseif strcmpi(corrMethodStr,'Kendall')
        corrMethod = 'Kendall';
    end
    [rho,pval] = corr(xMat,yMat,'type',corrMethod);
    [rho_pval] = [diag(rho),diag(pval)];
    n = size(xMat,1);
    [critRho10,critRho05,critRho01] = getCriticalTvalue(n);
end

function [critRho10,critRho05,critRho01] = getCriticalTvalue(sampleSizeArray)
    % use tinv to compute crtical tval for 0.1,0.05, and 0.01 pval
    % compute the critical rho vals for the pVals = 0.1,0.05,0.01 to test
    % for significance
    % use t = r*sqrt((n-2)/(1-r^2)) for t value
    n = sampleSizeArray;
    % these are tied to var names rho10,rho05,rho01
    levels = [0.1,0.05,0.01];
    tCrit = arrayfun(@(x) tinv(levels,x),n,'UniformOutput',false);
    rhoCrit = arrayfun(@(x) sqrt((tCrit{x}.^2)./(n(x)-2+tCrit{x}.^2)),(1:numel(tCrit))','UniformOutput',false);
    [critRho10,critRho05,critRho01] =cellfun(@(x) deal(x(1),x(2),x(3)),rhoCrit,'UniformOutput',false);
    critRho10 = cell2mat(critRho10);
    critRho05 = cell2mat(critRho05);
    critRho01 = cell2mat(critRho01);
end
