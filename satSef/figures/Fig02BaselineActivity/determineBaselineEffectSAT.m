function [ varargout ] = determineBaselineEffectSAT( binfo , ninfo , nstats , spikes , varargin )
%determineBaselineEffectSAT This function tests for a significant effect of
%SAT condition on baseline discharge rate for single neurons. That actual
%test is the (non-parametric) Mann-Whitney U-test.
% 

args = getopt(varargin, {{'area=','SEF'}, {'monkey=',{'D','E','Q','S'}}});

idxArea = ismember({ninfo.area}, args.area);
idxMonkey = ismember({ninfo.monkey}, args.monkey);

idxVis = ([ninfo.visGrade] >= 2);   idxMove = ([ninfo.moveGrade] >= 2);
idxErr = ([ninfo.errGrade] >= 2);   idxRew = (abs([ninfo.rewGrade]) >= 2);

idxKeep = (idxArea & idxMonkey & (idxVis | idxMove | idxErr | idxRew));

NUM_CELLS = sum(idxKeep);
ninfo = ninfo(idxKeep);
spikes = spikes(idxKeep);

T_BASE  = 3500 + [-600, 20];

for cc = 1:NUM_CELLS
  kk = ismember({binfo.session}, ninfo(cc).sess);
  ccNS = ninfo(cc).unitNum; %index nstats correctly
  
  %index by isolation quality
  idxIso = identify_trials_poor_isolation_SAT(ninfo(cc), binfo(kk).num_trials, 'task','SAT');
  %index by trial outcome
  idxCorr = ~(binfo(kk).err_dir | binfo(kk).err_time | binfo(kk).err_nosacc | binfo(kk).err_hold);
  %index by condition
  trialAcc = find((binfo(kk).condition == 1) & idxCorr & ~idxIso);
  trialFast = find((binfo(kk).condition == 3) & idxCorr & ~idxIso);
  
  nTrialAcc = length(trialAcc);
  nTrialFast = length(trialFast);
  
  spkCtAcc = NaN(1,nTrialAcc);
  for jj = 1:nTrialAcc
    spkTime_jj = spikes(cc).SAT{trialAcc(jj)};
    spkCtAcc(jj) = sum((spkTime_jj > T_BASE(1)) & (spkTime_jj < T_BASE(2)));
  end%for:trialAccurate(jj)
  
  spkCtFast = NaN(1,nTrialFast);
  for jj = 1:nTrialFast
    spkTime_jj = spikes(cc).SAT{trialFast(jj)};
    spkCtFast(jj) = sum((spkTime_jj > T_BASE(1)) & (spkTime_jj < T_BASE(2)));
  end%for:trialFast(jj)
  
  %Mann-Whitney U test for the difference between conditions (independent samples)
  [~,hSig,tmp] = ranksum(spkCtFast, spkCtAcc, 'alpha',0.05);
  if (hSig == 1)
    if (tmp.zval < 0) %Acc > Fast
      nstats(ccNS).blineEffect = -1;
    else %Fast > Acc
      nstats(ccNS).blineEffect = 1;
    end
  else%no SAT effect on baseline
    nstats(ccNS).blineEffect = 0;
  end
end%for:cells(cc)

if (nargout > 0)
  varargout{1} = nstats;
end

%% Output
nstats = nstats(idxKeep);
nFgA = sum(([nstats.blineEffect] ==  1));
nAgF = sum(([nstats.blineEffect] == -1));

fprintf('Number of neurons with Fast > Acc: %d/%d\n', nFgA, NUM_CELLS)
fprintf('Number of neurons with Acc > Fast: %d/%d\n', nAgF, NUM_CELLS)

end%fxn:determineBaselineEffectSAT()
