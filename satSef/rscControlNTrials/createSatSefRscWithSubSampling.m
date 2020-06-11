function [spkCorr] = createSatSefRscWithSubSampling()
% for figure8: use: satSef/figures/Fig08SpikeCorr/createSpikeCorrWithSubSampling.m
createSatSefDb = 'satSef/figures/Fig08SpikeCorr/createSpikeCorrWithSubSampling.m';
if (1)
    
    warning('for figure8 of SAT-SEF paper use: %s',createSatSefDb);
    edit('satSef/figures/Fig08SpikeCorr/createSpikeCorrWithSubSampling.m')
    error('Open createSatSefRscWithSubSamplingV2.m and read comments!')
end
% Create spike count correlation data set
% for pairs of units recorded from same session for cross areas.
% FOR-each-session in sessions DO --> cant do because for some units in the
% session a varying no. of trials have to be removed, bu tnot for other
% units. So we do not have a 'clear' nTrialsPerOutcome that is same for all
% recorded units.  Instead we need to do the computations on a
% neuron-neuron-pair basis. This is good because we do not have to
% eliminate whole session.
% crossPairs = getAllCrossPairs-no-filtering
% FOR-each-crossPair in crossPairs DO 
%   trialNos4SatConds = getTrialNosForAllSatConds(removeTrialNos-
%                             for-pair-if-present-on-any-unit-in-pair)
%   nTrials4SatConds = countTrials4SatConds(trialNos4SatConds)
%   IF min(trialNos4SatConds) < nTrialsThreshold THEN
%       ignore crossPair, go to next crossPair
%   ELSE-IF-min(trialNos4SatConds) >= nTrialsThreshold THEN
%       nTrials4SubSample = min(nTrials4SatConds)
%       satCondWithMinTrials = whichSatCondIs(nTrials4SubSample)
%       satConds2SubSample = whichSatCondIs(>nTrials4SubSample)
%       FOR-each-satCond in allSatConds DO
%           trialsNos4SatCond = getFrom(trialNos4SatConds,satCond)
%           nTrials4SatCond = count(trialsNos4SatCond)
%           [X_unitSpkCountByTrial, Y_unitSpkCountByTrial] =
%                getSpkCounts(X_unitAlignedTimes,Y_unitAlignedTimes,timWin)
%           [X_unitSpkCountMean,X_unitSpkCountStd] =
%                getSpkCountStats(X_unitAlignedSpikes)
%           [Y_unitSpkCountMean,Y_unitSpkCountStd] =
%                getSpkCountStats(Y_unitAlignedSpikes)
%           [cpRsc,cpPairPval] = corr(X_unitSpkCountByTrial,
%                Y_unitSpkCountByTrial) Pearson's
%            IF-is-member(satConds2SubSample,satCond) THEN
%                FOR-1-to-nSubSamples DO
%                    subSampIdx = subsample(1:nTrials4SatCond,
%                        nTrials4SubSample) 
%                    estimatedRscVec(loop-count) =
%                        corr(X_unitSpkCountByTrial(subSampIdx), 
%                        Y_unitSpkCountByTrial(subSampIdx)) Pearson's
%                LOOP-FOR-nSubSamples
%                estimatedRsc4Cond = mean(estimatedRscVec)
%                confIntRsc4Cond = CI(estimatedRscVec,0.95) 
%            ELSE-IF-not-is-member(satConds2SubSample,satCond) THEN
%                GOTO-next-satCond
%            END-IF- is-member(satConds2SubSample,satCond)
%       LOOP-FOR-eachSatCond
%   END-ELSE-IF-min(nTrial count) >= nTrialsThreshold
% LOOP-FOR-each-crossPair

% Expects the following files in specified location:
% 1. Inofrmation about all possible cell pairs:
%        'dataProcessed/satSefPaper/dataset/SAT_SEF_PAIR_CellInfoDB.mat'
% 2. Trial types of all sessions (Accurate, Fast, Correct,...):
%        'dataProcessed/satSefPaper/dataset/SAT_SEF_TrialTypesDB.mat'
% 3. Event times for all trials and all sessions:
%         'dataProcessed/satSefPaper/dataset/SAT_SEF_TrialEventTimesDB.mat'
% 4. Get resptime and set it as SaccadePrimaryTempo
%         'dataProcessed/dataset/binfo_moves_SAT.mat'
% 5. Spike time data for all units of all sessions:
%         'dataProcessed/dataset/spikes_SAT.mat'
%

% fx_subsample = @(vec,n) datasample(vec,n,'Replace',false);

%% File refs for data to computing Rsc
    warning('off');
    %Options for Spk Corr computation
    datasetDir = 'dataProcessed/satSefPaper/dataset';
    % Files for getting data to compute Rsc
    pairsFile = fullfile(datasetDir,'SAT_SEF_PAIR_CellInfoDB.mat');
    trialTypesFile = fullfile(datasetDir,'SAT_SEF_TrialTypesDB.mat');
    trialEventTimesFile = fullfile(datasetDir,'SAT_SEF_TrialEventTimesDB.mat');
    spikeTimesFile = fullfile(datasetDir,'spikes_SAT.mat');
    % output file
    outFile = 'newRscWithSubSampling.mat';

    % alignment
    % Setup time windows for different event time alignment, the field names
    % SHALL correspond to column names for trialEventTimes below.
    alignNames = {'PostSaccade'};
    alignEvents = {'SaccadePrimary'};
    alignTimeWin = {[-100 500]};

    conditions = {'AccurateCorrect';'AccurateErrorChoice';'AccurateErrorTiming';
        'FastCorrect';    'FastErrorChoice';    'FastErrorTiming'
        };
    staticWinsAlignTs(1).PostSaccade = [0 150];
    % minimum number of trials for all conditions, if not, then ignore pair  
    nTrialsThreshold = 10; %10;%3 should reproduce the old result 
    nSubSamples = 1000; % number of times to subsample
    
    %% Process  (use parfor, when available)
    p = gcp('nocreate');
    nThreads = 0;
    if ~isempty(p)
        nThreads = p.NumWorkers;
    end
    
    %%  Load data needed to compute Rsc
    crossPairs = getCrossAreaPairs(pairsFile);
    spikesSat = getSpikes(spikeTimesFile);
    sessionTrialTypes = getSessionTrialTypes(trialTypesFile);
    sessionEventTimes = getSessionEventTimes(trialEventTimesFile);

    %%
    tic
    nCrossPairs = size(crossPairs,1);
    spkCorr = table();
    pctRunOnAll warning off;
    %parfor (cp = 1:nCrossPairs,nThreads)%nCrossPairs
    parfor (cp = 1:nCrossPairs,nThreads)%nCrossPairs
        %%
        opts = struct();
        crossPair = crossPairs(cp,:);
        sess = crossPair.X_sess{1};
        trialTypes = sessionTrialTypes(ismember(sessionTrialTypes.session,sess),:); %#ok<*PFBNS>
        trRem = getTrialNosToRemove(crossPair);
        [trialNos4SatConds, nTrials4SatConds] = ...
            getTrialNosForAllSatConds(trialTypes,trRem,conditions);
        nTrials4SubSample = min(struct2array(nTrials4SatConds));
        if nTrials4SubSample < nTrialsThreshold
            % ignore this pair and go to next pair
            continue;
        end
        % spike times
        xSpkTimes = spikesSat{crossPair.X_unitNum}';
        ySpkTimes = spikesSat{crossPair.Y_unitNum}';
        evntTimes = sessionEventTimes(ismember(sessionEventTimes.session,sess),:);
        for sc = 1:numel(conditions)
            condition = conditions{sc};
            selTrials = trialNos4SatConds.(condition);
            selTrialsSorted = selTrials; % no sorting
            for evId = 1:numel(alignEvents)
                alignedName = alignNames{evId};
                alignedEvent = alignEvents{evId};
                alignedTimeWin = alignTimeWin{evId};
                alignTime = evntTimes.CueOn{1};
                if ~strcmp(alignedEvent,'CueOn')
                    alignTime = alignTime + double(evntTimes.(alignedEvent){1}(:));
                end  
                alignTime = alignTime(selTrialsSorted);
                XAligned = SpikeUtils.alignSpikeTimes(xSpkTimes(selTrialsSorted),alignTime, alignedTimeWin);
                YAligned = SpikeUtils.alignSpikeTimes(ySpkTimes(selTrialsSorted),alignTime, alignedTimeWin);
                tempRast = SpikeUtils.rasters(XAligned,alignedTimeWin);
                XRasters = tempRast.rasters;
                tempRast = SpikeUtils.rasters(YAligned,alignedTimeWin);
                YRasters = tempRast.rasters;
                rasterBins = tempRast.rasterBins;                               
                
                staticWin = staticWinsAlignTs(1).(alignedName);
                fieldSuffix = num2str(range(staticWin),'_%dms'); 
                xSpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{XRasters},{staticWin},'UniformOutput',false);
                xMeanFrWin = mean(xSpkCounts{1})*1000/range(staticWin);
                ySpkCounts = cellfun(@(r,x,w) sum(x(:,r>=w(1) & r<=w(2)),2),...
                    {rasterBins},{YRasters},{staticWin},'UniformOutput',false);
                yMeanFrWin = mean(ySpkCounts{1})*1000/range(staticWin);
                [rho_pval] = getSpikeCountCorr(xSpkCounts{1},ySpkCounts{1},'Pearson');
                % Do sub-sampling to estimate Rho and CI
                % nTrials2SubSample, selTrials
                [rhoEst,rhoEstSem,percentileCI,normalCI] = ...
                     getEstimatedRhoAndConfInterval(xSpkCounts{1},...
                                                    ySpkCounts{1},...
                                                    nTrials4SubSample,...
                                                    nSubSamples);

                % output table/struct
                % add crosspair info
                cN = getPairColNmes();
                cpTemp = table2struct(crossPair,'ToScalar',true);
                for c = 1:numel(cN)
                    opts(evId,sc).(cN{c}) = cpTemp.(cN{c});
                end
                opts(evId,sc).condition = condition;
                opts(evId,sc).alignedName = alignedName;
                opts(evId,sc).alignedEvent = alignedEvent;
                opts(evId,sc).alignedTimeWin = {alignedTimeWin};
                opts(evId,sc).alignTime = alignTime;
                opts(evId,sc).trialNosByCondition = selTrialsSorted;
                opts(evId,sc).nTrials = numel(selTrialsSorted);

                opts(evId,sc).(['xSpkCount_win' fieldSuffix]) = xSpkCounts;
                opts(evId,sc).(['ySpkCount_win' fieldSuffix]) = ySpkCounts;
                opts(evId,sc).(['xMeanFr_spkPerSec_win' fieldSuffix]) = xMeanFrWin;
                opts(evId,sc).(['yMeanFr_spkPerSec_win' fieldSuffix]) = yMeanFrWin;

                opts(evId,sc).(['rho_pval_win' fieldSuffix]) = {staticWin};                    
                opts(evId,sc).(['rhoRaw' fieldSuffix]) = rho_pval(1);
                opts(evId,sc).(['pvalRaw' fieldSuffix]) = rho_pval(2);
                opts(evId,sc).(['signifRaw_05' fieldSuffix]) = rho_pval(2) < 0.05;  

                opts(evId,sc).nTrials4SubSample = nTrials4SubSample;
                opts(evId,sc).nSubSamples = nSubSamples;
                opts(evId,sc).(['rhoEstRaw' fieldSuffix]) = rhoEst;
                opts(evId,sc).(['rhoEstSem' fieldSuffix]) = rhoEstSem;
                opts(evId,sc).(['normalCI95' fieldSuffix]) = {normalCI};
                opts(evId,sc).(['percentileCI95' fieldSuffix]) = {percentileCI};
                opts(evId,sc).rhoRawWithinNormalCI = normalCI(1) < rho_pval(1) & rho_pval(1) < normalCI(2);
                opts(evId,sc).rhoRawWithinPercentileCI = percentileCI(1) < rho_pval(1) & rho_pval(1) < percentileCI(2);               
                
            end % for aligned event
            %spkCorr = [spkCorr; [crossPair(1,getPairColNmes) struct2table(opts,'AsArray',true)]];
        end % for condition
        spkCorr = [spkCorr;struct2table(opts)];
        fprintf('Done pair %d of %d\n',cp,nCrossPairs)
    end
toc
save(outFile,'-v7.3','spkCorr');
end

function [rhoEst,rhoEstSem,percentileCI,normalCI] = getEstimatedRhoAndConfInterval(xSpkCounts,ySpkCounts,nTrials4SubSample,nSubSamples)
    % inline fx for subsampling see DATASAMPLE
    subSampleIdxs = arrayfun(@(x) datasample(1:numel(xSpkCounts),nTrials4SubSample,'Replace',true)',(1:nSubSamples),'UniformOutput',false);
    temp = cellfun(@(x) getSpikeCountCorr(xSpkCounts(x),ySpkCounts(x),'Pearson'),subSampleIdxs','UniformOutput',false);
    temp = cell2mat(temp);
    rhoVec = temp(:,1);
    % compute mean & sem
    rhoEst = mean(rhoVec);
    rhoEstSem = std(rhoVec)/sqrt(numel(rhoVec));
    % compute t-statistic for 0.025, 0.975
    ts = tinv([0.025,0.975],nTrials4SubSample-1);
    normalCI = rhoEst + rhoEstSem*ts;
    percentileCI =  prctile(rhoVec,[2.5, 97.5]);
end


function [cellPairs] = getCrossAreaPairs(pairsFile)
    allCellPairs = load(pairsFile);
    allCellPairs = allCellPairs.satSefPairCellInfoDB;
     monkIdsToDo = {'D','E'};
    allCellPairs = allCellPairs(ismember([allCellPairs.X_monkey],monkIdsToDo),:);
    idxCrossArea = ismember(allCellPairs.X_area,'SEF') & ...
                   (ismember(allCellPairs.Y_area,'FEF') | ...
                    ismember(allCellPairs.Y_area,'SC'));
    cellPairs = allCellPairs(idxCrossArea,:);
    assert(isequal(cellPairs.X_sess,cellPairs.Y_sess),'********Fatal: Error X-Unit sessions and Y-Unit sessions do not match');
    fprintf('Done getCrossAreaPairs()\n')
end

function [spikesSat] = getSpikes(spikeTimesFile)
    spikesSat = load(spikeTimesFile);
    spikesSat = spikesSat.spikesSAT;
    fprintf('Done getSpikes()\n')
end

function [sessionTrialTypes] = getSessionTrialTypes(trialTypesFile)
    sessionTrialTypes = load(trialTypesFile);
    sessionTrialTypes = sessionTrialTypes.TrialTypesDB;
    fprintf('Done getSessionTrialTypes()\n')
end

function [sessionEventTimes] = getSessionEventTimes(trialEventTimesFile)
    sessionEventTimes = load(trialEventTimesFile);
    sessionEventTimes = sessionEventTimes.TrialEventTimesDB;
    fprintf('Done getSessionEventTimes()\n')
end

function [trRem] = getTrialNosToRemove(crossPair)
    trRem = [crossPair.X_trRemSAT{1};crossPair.Y_trRemSAT{1}];
    if ~isempty(trRem)
        temp = [];
        for ii = 1:size(trRem,1)
            temp = [temp [trRem(ii,1):trRem(ii,2)]]; %#ok<NBRAK,AGROW>
        end
        trRem = unique(temp(:));
    end
    %fprintf('Done getTrialNosToRemove()\n')
end

function [trialNos4SatConds,nTrials4SatConds] = getTrialNosForAllSatConds(trialTypes,trRem,conditions)
   trialNos4SatConds = struct();
   nTrials4SatConds = struct();
   for c = 1:numel(conditions)
       condition = conditions{c};
       temp = trialTypes.(condition){1};
       temp(trRem) = 0;
       temp = find(temp);       
       trialNos4SatConds.(condition) = temp;
   end
   % remove trials common to *ErrorChoice and *ErrorTiming
   % May not need this since they already seem to be mutually exclusive
   commonTrialNos = intersect(trialNos4SatConds.AccurateErrorChoice,trialNos4SatConds.AccurateErrorTiming);
   trialNos4SatConds.AccurateErrorChoice = setdiff(trialNos4SatConds.AccurateErrorChoice,commonTrialNos);
   trialNos4SatConds.AccurateErrorTiming = setdiff(trialNos4SatConds.AccurateErrorTiming,commonTrialNos);
   commonTrialNos = intersect(trialNos4SatConds.FastErrorChoice,trialNos4SatConds.FastErrorTiming);
   trialNos4SatConds.FastErrorChoice = setdiff(trialNos4SatConds.FastErrorChoice,commonTrialNos);
   trialNos4SatConds.FastErrorTiming = setdiff(trialNos4SatConds.FastErrorTiming,commonTrialNos);
   % nTrialsFor sat conds
   for c = 1:numel(conditions)
       condition = conditions{c};
       nTrials4SatConds.(condition) = numel(trialNos4SatConds.(condition));
   end   
   %fprintf('Done getTrialNosForAllSatConds()\n')   
end

function [colNames] = getPairColNmes()
    % cellPairInfo fields to extract
    colNames = {
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
        };
end