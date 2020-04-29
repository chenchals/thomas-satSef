function [ ] = Fig01B_Plot_ErrRate_X_RT( binfo , pSacc )
%Fig01B_Plot_ErrRate_X_RT() Summary of this function goes here
%   Detailed explanation goes here

PLOT = true;
STATS = true;

%isolate sessions from MONKEY
MONKEY = {'D','E'};         sessKeep = ismember(binfo.monkey, MONKEY);
NUM_SESS = sum(sessKeep);   binfo = binfo(sessKeep, :);   pSacc = pSacc(sessKeep, :);


%% Compute RT/ER and split on Task Condition
errRate_Acc = NaN(1,NUM_SESS);   rt_Acc = NaN(1,NUM_SESS);
errRate_Fast = NaN(1,NUM_SESS);  rt_Fast = NaN(1,NUM_SESS);

for kk = 1:NUM_SESS
  
  %index by condition
  idxAcc = (binfo.condition{kk} == 1) & ~isnan(binfo.deadline{kk});
  idxFast = (binfo.condition{kk} == 3) & ~isnan(binfo.deadline{kk});
  %index by trial outcome
  idxCorr = ~(binfo.err_dir{kk} | binfo.err_time{kk} | binfo.err_nosacc{kk});
  idxErr = (binfo.err_dir{kk});
  
  rt_Acc(kk) = median(pSacc.resptime{kk}(idxAcc & idxCorr));
  rt_Fast(kk) = median(pSacc.resptime{kk}(idxFast & idxCorr));
  errRate_Acc(kk) = sum(idxAcc & idxErr) / sum(idxAcc);
  errRate_Fast(kk) = sum(idxFast & idxErr) / sum(idxFast);
  
end%for:session(kk)


%% Split RT/ER on Search Difficulty
idxMore = (binfo.taskType == 1); NUM_MORE = sum(idxMore); %more efficient
idxLess = (binfo.taskType == 2); NUM_LESS = sum(idxLess); %less efficient

%split RT by condition and efficiency
rt_AccMore = rt_Acc(idxMore);         rt_AccLess = rt_Acc(idxLess);
er_AccMore = errRate_Acc(idxMore);    er_AccLess = errRate_Acc(idxLess);
rt_FastMore = rt_Fast(idxMore);       rt_FastLess = rt_Fast(idxLess);
er_FastMore = errRate_Fast(idxMore);  er_FastLess = errRate_Fast(idxLess);

%compute mean and SE of response time
mu.RT_AM = mean(rt_AccMore);    se.RT_AM = std(rt_AccMore)/sqrt(NUM_MORE);
mu.RT_AL = mean(rt_AccLess);    se.RT_AL = std(rt_AccLess)/sqrt(NUM_LESS);
mu.RT_FM = mean(rt_FastMore);   se.RT_FM = std(rt_FastMore)/sqrt(NUM_MORE);
mu.RT_FL = mean(rt_FastLess);   se.RT_FL = std(rt_FastLess)/sqrt(NUM_LESS);
%compute mean and SE of error rate
mu.ER_AM = mean(er_AccMore);    se.ER_AM = std(er_AccMore)/sqrt(NUM_MORE);
mu.ER_AL = mean(er_AccLess);    se.ER_AL = std(er_AccLess)/sqrt(NUM_LESS);
mu.ER_FM = mean(er_FastMore);   se.ER_FM = std(er_FastMore)/sqrt(NUM_MORE);
mu.ER_FL = mean(er_FastLess);   se.ER_FL = std(er_FastLess)/sqrt(NUM_LESS);


%% Plotting
if (PLOT)
  figure()

  subplot(1,3,1); hold on
  errorbarxy([mu.RT_FM mu.RT_AM], [mu.ER_FM mu.ER_AM], [se.RT_FM se.RT_AM], [se.ER_FM se.ER_AM], {'k-','k','k'})
  errorbarxy([mu.RT_FL mu.RT_AL], [mu.ER_FL mu.ER_AL], [se.RT_FL se.RT_AL], [se.ER_FL se.ER_AL], {'k-','k','k'})
  xlim([245 600]); ylim([.05 .45])

  subplot(1,3,2); hold on
  errorbar([mu.RT_FM mu.RT_AM], [se.RT_FM se.RT_AM], 'CapSize',0, 'LineWidth',1, 'Color','k')
  errorbar([mu.RT_FL mu.RT_AL], [se.RT_FL se.RT_AL], 'CapSize',0, 'LineWidth',2, 'Color','k')
  xlim([0.9 2.1]); xticks([1 2]); xticklabels({'Fast','Accurate'})

  subplot(1,3,3); hold on
  errorbar([mu.ER_FM mu.ER_AM], [se.ER_FM se.ER_AM], 'CapSize',0, 'LineWidth',1, 'Color','k')
  errorbar([mu.ER_FL mu.ER_AL], [se.ER_FL se.ER_AL], 'CapSize',0, 'LineWidth',2, 'Color','k')
  xlim([0.9 2.1]); xticks([1 2]); xticklabels({'Fast','Accurate'})

  ppretty([8,1.8])
end

%% Stats -- Two-way between-subjects ANOVA
if (STATS)
    RT = [rt_AccMore rt_AccLess rt_FastMore rt_FastLess]';
    ER = [er_AccMore er_AccLess er_FastMore er_FastLess]';
    F_Condition = [ones(1,NUM_SESS) 2*ones(1,NUM_SESS)]';
    F_Efficiency = [ones(1,NUM_MORE) 2*ones(1,NUM_LESS) ones(1,NUM_MORE) 2*ones(1,NUM_LESS)]';
    F_Session = [(1:NUM_SESS) (1:NUM_SESS)];
    
    % try new method for all ANOVA
    anovaTbl = table();
    % Value
    anovaTbl.RT = RT;
    % 2-way anova 
    % Factor-1 - categorical use names
    anovaTbl.Condition = arrayfun(@(x) sprintf('Condition_%d',x), F_Condition,'UniformOutput',false);
    % Factor-2 - use names
    anovaTbl.Efficiency =  arrayfun(@(x) sprintf('Efficiency%d',x), F_Efficiency,'UniformOutput',false);
    % do anova
    anovaModel = 'full'; %  [linear|interaction|full]
    multiCompareFlag = 1; % [0|1]
    alpha = 0.05;
    [anovaResults] = satAnova(anovaTbl,'model','full','doMultiCompare',true,'alpha',0.05)
end

end % fxn :: Fig01B_Plot_ErrRate_X_RT()


