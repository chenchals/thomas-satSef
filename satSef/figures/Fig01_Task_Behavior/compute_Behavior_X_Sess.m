function [ ] = compute_Behavior_X_Sess( binfo , pSacc , sSacc )
%compute_Behavior_X_Sess Summary of this function goes here
%   Detailed explanation goes here

%isolate sessions from MONKEY
MONKEY = {'D','E'};         sessKeep = ismember(binfo.monkey, MONKEY);
NUM_SESS = sum(sessKeep);   binfo = binfo(sessKeep, :);   pSacc = pSacc(sessKeep, :);   sSacc = sSacc(sessKeep, :);

%% Initializations

dlineAcc = NaN(1,NUM_SESS);    RTAcc = NaN(1,NUM_SESS);
dlineFast = NaN(1,NUM_SESS);   RTFast = NaN(1,NUM_SESS);

PerrChcAcc = NaN(1,NUM_SESS);  PerrTimeAcc = NaN(1,NUM_SESS);
PerrChcFast = NaN(1,NUM_SESS); PerrTimeFast = NaN(1,NUM_SESS);

isiAcc = NaN(1,NUM_SESS);
isiFast = NaN(1,NUM_SESS);

%% Collect data

for kk = 1:NUM_SESS
  
  idxAcc = (binfo.condition{kk} == 1) & ~isnan(binfo.deadline{kk});
  idxFast = (binfo.condition{kk} == 3) & ~isnan(binfo.deadline{kk});
  
  idxCorr = ~(binfo.err_dir{kk} | binfo.err_time{kk} | binfo.err_nosacc{kk});
  idxErrChc = binfo.err_dir{kk} & ~binfo.err_time{kk};
  idxErrTime = binfo.err_time{kk} & ~binfo.err_dir{kk};
  
  dlineAcc(kk) = median(binfo.deadline{kk}(idxAcc));
  dlineFast(kk) = median(binfo.deadline{kk}(idxFast));
  
  RTAcc(kk) = median(pSacc.resptime{kk}(idxAcc & idxCorr));
  RTFast(kk) = median(pSacc.resptime{kk}(idxFast & idxCorr));
  
  PerrChcAcc(kk) = sum(idxAcc & idxErrChc) / sum(idxAcc);
  PerrChcFast(kk) = sum(idxFast & idxErrChc) / sum(idxFast);
  
  PerrTimeAcc(kk) = sum(idxAcc & idxErrTime) / sum(idxAcc);
  PerrTimeFast(kk) = sum(idxFast & idxErrTime) / sum(idxFast);
  
  ISIkk = double(sSacc.resptime{kk}) - (double(pSacc.resptime{kk}) + double(pSacc.duration{kk}));
  idxNoPP = (sSacc.resptime{kk} == 0);
  
  isiAcc(kk) = median(ISIkk(idxAcc & idxErrChc & ~idxNoPP));
  isiFast(kk) = median(ISIkk(idxFast & idxErrChc & ~idxNoPP));
  
end%for:session(kk)

%% Print mean +/- SE
fprintf('Response deadline Acc: %g +/- %g\n', mean(dlineAcc), std(dlineAcc)/sqrt(NUM_SESS))
fprintf('Response deadline Fast: %g +/- %g\n\n', mean(dlineFast), std(dlineFast)/sqrt(NUM_SESS))

fprintf('RT Acc: %g +/- %g\n', mean(RTAcc), std(RTAcc)/sqrt(NUM_SESS))
fprintf('RT Fast: %g +/- %g\n', mean(RTFast), std(RTFast)/sqrt(NUM_SESS))

fprintf('ISI Acc: %g +/- %g\n', mean(isiAcc), std(isiAcc)/sqrt(NUM_SESS))
fprintf('ISI Fast: %g +/- %g\n', mean(isiFast), std(isiFast)/sqrt(NUM_SESS))

fprintf('P[ChcErr] Acc: %g +/- %g\n', mean(PerrChcAcc), std(PerrChcAcc)/sqrt(NUM_SESS))
fprintf('P[ChcErr] Fast: %g +/- %g\n', mean(PerrChcFast), std(PerrChcFast)/sqrt(NUM_SESS))

fprintf('P[TimeErr] Acc: %g +/- %g\n', mean(PerrTimeAcc), std(PerrTimeAcc)/sqrt(NUM_SESS))
fprintf('P[TimeErr] Fast: %g +/- %g\n', mean(PerrTimeFast), std(PerrTimeFast)/sqrt(NUM_SESS))

fprintf('\n')

end%fxn:compute_Behavior_X_Sess()
