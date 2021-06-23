function [ ] = computeRankRho_BaselineCount_X_RT_SAT( binfo , moves , ninfo , spikes , varargin )
%computeRankRho_BaselineCount_X_RT_SAT Summary of this function goes here
%   Detailed explanation goes here

args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E','Q','S'}}});

if ~any(ismember(args.monkey, {'Q','S'}))
  binfo = binfo(1:16);
end

idxArea = ismember({ninfo.area}, args.area);
idxMonkey = ismember({ninfo.monkey}, args.monkey);

idxVis = ([ninfo.visGrade] >= 2);   idxMove = ([ninfo.moveGrade] >= 2);
idxErr = ([ninfo.errGrade] >= 2);   idxRew = (abs([ninfo.rewGrade]) >= 2);
idxTaskRel = (idxVis | idxMove | idxErr | idxRew);

idxKeep = (idxArea & idxMonkey & idxTaskRel);

NUM_CELLS = sum(idxKeep);
ninfo = ninfo(idxKeep);
spikes = spikes(idxKeep);

RTLIM_ACC = [390 800];
RTLIM_FAST = [150 450];

T_BLINE = 3500 + [-600 20];

%initializations
rhoSpearmanAcc = NaN(1,NUM_CELLS);    pvalSpearmanAcc = NaN(1,NUM_CELLS);
rhoSpearmanFast = NaN(1,NUM_CELLS);   pvalSpearmanFast = NaN(1,NUM_CELLS);

for cc = 1:NUM_CELLS
  fprintf('Unit %s - %s\n', ninfo(cc).sess, ninfo(cc).unit);
  
  %compute spike count for all trials
  spkCtCC = cellfun(@(x) sum((x > T_BLINE(1)) & (x < T_BLINE(2))), spikes(cc).SAT);
  
  kk = ismember({binfo.session}, ninfo(cc).sess);
  RTkk = double(moves(kk).resptime);
  
  %index by isolation quality
  idxIso = identify_trials_poor_isolation_SAT(ninfo(cc), binfo(kk).num_trials);
  %index by trial outcome
  idxCorr = ~(binfo(kk).err_dir | binfo(kk).err_time | binfo(kk).err_nosacc | binfo(kk).err_hold);
  %index by condition and RT limits
  idxAcc = ((binfo(kk).condition == 1) & idxCorr & ~idxIso & ~(RTkk < RTLIM_ACC(1) | RTkk > RTLIM_ACC(2) | isnan(RTkk)));
  idxFast = ((binfo(kk).condition == 3) & idxCorr & ~idxIso & ~(RTkk < RTLIM_FAST(1) | RTkk > RTLIM_FAST(2) | isnan(RTkk)));
  
  spkCountAcc = spkCtCC(idxAcc);    RTacc = RTkk(idxAcc);
  spkCountFast = spkCtCC(idxFast);  RTfast = RTkk(idxFast);
  
  %compute Spearman rank correlation coefficient
  [rAcc,pAcc] = corr(RTacc', spkCountAcc', 'Type','Spearman');
  [rFast,pFast] = corr(RTfast', spkCountFast', 'Type','Spearman');
  
  rhoSpearmanAcc(cc) = rAcc;      pvalSpearmanAcc(cc) = pAcc;
  rhoSpearmanFast(cc) = rFast;    pvalSpearmanFast(cc) = pFast;
  
end%for:cell(cc)

%separate neurons by level of search efficiency
idxMore = [ninfo.taskType] == 1;  NUM_MORE = sum(idxMore);
idxLess = [ninfo.taskType] == 2;  NUM_LESS = sum(idxLess);

%find units with a significant correlation
ALPHA = 0.06;
idxAccPos = (pvalSpearmanAcc <= ALPHA & rhoSpearmanAcc > 0);
idxAccNeg = (pvalSpearmanAcc <= ALPHA & rhoSpearmanAcc < 0);
idxFastPos = (pvalSpearmanFast <= ALPHA & rhoSpearmanFast > 0);
idxFastNeg = (pvalSpearmanFast <= ALPHA & rhoSpearmanFast < 0);

fprintf('More efficient:\n')
fprintf('Acc: (+) = %d  (-) = %d  / %d\n', sum(idxAccPos & idxMore), sum(idxAccNeg & idxMore), NUM_MORE)
fprintf('Fast: (+) = %d  (-) = %d  / %d\n', sum(idxFastPos & idxMore), sum(idxFastNeg & idxMore), NUM_MORE)
fprintf('Less efficient:\n')
fprintf('Acc: (+) = %d  (-) = %d  / %d\n', sum(idxAccPos & idxLess), sum(idxAccNeg & idxLess), NUM_LESS)
fprintf('Fast: (+) = %d  (-) = %d  / %d\n', sum(idxFastPos & idxLess), sum(idxFastNeg & idxLess), NUM_LESS)

end%fxn:computeRankRho_BaselineCount_X_RT_SAT()
