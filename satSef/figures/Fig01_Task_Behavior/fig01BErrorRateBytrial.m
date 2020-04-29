function [] = fig01BErrorRateBytrial(binfo, pSacc)
%FIG01BErrorRateBYTRIAL Summary of this function goes here
%   Detailed explanation goes here
PLOT = true;
STATS = true;

%% Use sessions from monkeys
monkeys = {'D','E'}; 
idxMonks = ismember(binfo.monkey,monkeys);   
binfo = binfo(idxMonks,:);   
pSacc = pSacc(idxMonks,:);
nSess = sum(idxMonks);

%% Create a table for Reation time, error rate, efficiency by session for SAT conditions
rtErTbl = table();
% binfo condition: 1= Accurate and 3 = Fast, 
% so use same sequence for sat conditions
satConditionVals = [1 3];
satConditions = {'Accurate','Fast'};
% idxMoreEfficient = (binfo.taskType == 1); %more efficient CS: less difficult
% idxLessEfficient = (binfo.taskType == 2); %less efficient CS: more difficult

for ss = 1:nSess
    sessNum = ss;
    % index for Corerect and error trials in session
    idxCorr = ~(binfo.err_dir{sessNum} | binfo.err_time{sessNum} | binfo.err_nosacc{sessNum});
    idxErr = binfo.err_dir{sessNum};
    if binfo.taskType(sessNum) == 1
        % idxMoreEfficient = (binfo.taskType == 1); %more efficient CS: less difficult
        efficiency ='MoreEfficient';
    elseif binfo.taskType(sessNum) == 2
        % idxLessEfficient = (binfo.taskType == 2); %less efficient CS: more difficult
        efficiency ='LessEfficient';    
    end    
    for con = 1:numel(satConditions)
        conditionVal = satConditionVals(con);
        satCondition = satConditions{con};
        idxSat = binfo.condition{ss} == conditionVal;
        % session
        temp = table();
        % temp.sessNum = sessNum;
        % Sat condition
        temp.satCondition = {satCondition};
        % Search efficiency
        temp.efficiency = {efficiency};
        % Reaction Time
        % rtAccurate(kk) = median(pSacc.resptime{kk}(idxAcc & idxCorr));
        temp.reactionTime =  median(pSacc.resptime{sessNum}(idxSat & idxCorr));
        % Error Rate
        % errRateAccurate(kk) = sum(idxAcc & idxErr) / sum(idxAcc);
        temp.errorRate = sum(idxSat & idxErr)/sum(idxSat);
        
        rtErTbl = [rtErTbl;temp];  %#ok<*AGROW>
    end
end
rtErTbl = sortrows(rtErTbl,{'efficiency','satCondition'},{'descend','descend'});

%% to check if OK
% the variable name given as first argument must exist in workspace
% check('rt_AccMore',rtErTbl([1:2:16],:))

%% compute mean and SEM for Reactions and Error Rate table and Plot results
rtErStatsTbl = grpstats(rtErTbl,{'satCondition','efficiency'},{'mean','sem'});
rtErStatsTbl = sortrows(rtErStatsTbl,{'efficiency','satCondition'},{'descend','descend'});
figure()
% Efficiency and errorRate
subplot(1,3,1);
hold on
idxM=ismember(rtErStatsTbl.efficiency,'MoreEfficient');
idxL=ismember(rtErStatsTbl.efficiency,'LessEfficient');
errorbar(rtErStatsTbl.mean_reactionTime(idxM),...
         rtErStatsTbl.mean_errorRate(idxM),...
         rtErStatsTbl.sem_errorRate(idxM),...
         '+-k','MarkerSize',8,'LineWidth',1,'CapSize',0)
errorbar(rtErStatsTbl.mean_reactionTime(idxL),...
         rtErStatsTbl.mean_errorRate(idxL),...
         rtErStatsTbl.sem_errorRate(idxL),...
         '+-k','MarkerSize',8,'LineWidth',1,'CapSize',0)
xlim([245 600]); ylim([.05 .40]);
xlabel('Reaction time (ms)')
ylabel('Error rate')


% condition and efficiency and reaction time  
subplot(1,3,2); 
hold on
errorbar([1 2],...
         rtErStatsTbl.mean_reactionTime(idxM),...
         rtErStatsTbl.sem_reactionTime(idxM),...
         '+-k','MarkerSize',1,'LineWidth',1,'CapSize',0)
errorbar([1 2],...
         rtErStatsTbl.mean_reactionTime(idxL),...
         rtErStatsTbl.sem_reactionTime(idxL),...
         '+-k','MarkerSize',1,'LineWidth',2,'CapSize',0)
 xlim([0.9 2.1]); xticks([1 2]); xticklabels({'Fast','Accurate'})
xlabel('SAT Condition')
ylabel('Reaction time (ms)')

% condition by efficiency and error rate 
subplot(1,3,3); 
hold on
errorbar([1 2],...
         rtErStatsTbl.mean_errorRate(idxM),...
         rtErStatsTbl.sem_errorRate(idxM),...
         '+-k','MarkerSize',1,'LineWidth',1,'CapSize',0)
errorbar([1 2],...
         rtErStatsTbl.mean_errorRate(idxL),...
         rtErStatsTbl.sem_errorRate(idxL),...
         '+-k','MarkerSize',1,'LineWidth',2,'CapSize',0)
xlim([0.9 2.1]); xticks([1 2]); xticklabels({'Fast','Accurate'})
xlabel('SAT Condition')
ylabel('Error rate')


end


%% Check values with previous code...
function [] = check(varname2Check,rtErTbl)
%%
% varname2Check must be in your workspace
var2Check = evalin('base',varname2Check);
% sat condition
if contains(varname2Check,'Acc')
    satCondition = 'Accurate';
elseif contains(varname2Check,'Fast')
    satCondition = 'Fast';
end
% efficiency
if endsWith(varname2Check,'More')
    efficiency = 'MoreEfficient';
elseif contains(varname2Check,'Fast')
    efficiency = 'LessEfficient';
end
% rt or er
if startsWith(varname2Check,'er_')
    rtOrEr = 'errorRate';
elseif startsWith(varname2Check,'rt_')
    rtOrEr = 'reactionTime';
end

zIdx=ismember(rtErTbl.satCondition,satCondition) & ismember(rtErTbl.efficiency,efficiency);
% assert(isequal(rtErTbl.(rtOrEr)(zIdx)',var2Check),sprintf('**%s does not agree!!\n',varname2Check));
if isequal(rtErTbl.(rtOrEr)(zIdx)',var2Check)
    fprintf('Verified OK\n');
else
    fprintf('*** Verification Failed for %s ***\n',varname2Check);
end

end
