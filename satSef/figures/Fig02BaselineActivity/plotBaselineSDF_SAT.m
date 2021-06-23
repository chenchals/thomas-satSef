function [ ] = plotBaselineSDF_SAT( binfo , moves , ninfo , nstats , spikes , varargin )
%plotBaselineSDF_SAT Summary of this function goes here
%   Detailed explanation goes here

args = getopt(varargin, {{'area=',{'SEF'}}, {'monkey=',{'D','E','Q','S'}}});

idxArea = ismember({ninfo.area}, args.area);
idxMonkey = ismember({ninfo.monkey}, args.monkey);

idxVis = ([ninfo.visGrade] >= 2);   idxMove = ([ninfo.moveGrade] >= 2);
idxErr = ([ninfo.errGrade] >= 2);   idxRew = (abs([ninfo.rewGrade]) >= 2);
idxTaskRel = (idxVis | idxMove);% | idxErr);% | idxRew);

idxKeep = (idxArea & idxMonkey & idxTaskRel);

NUM_CELLS = sum(idxKeep);
ninfo = ninfo(idxKeep);
nstats = nstats(idxKeep);
spikes = spikes(idxKeep);

T_STIM = 3500 + (-300 : 150); %from stimulus
N_SAMP = length(T_STIM); %number of samples consistent across epochs
T_FOCUSED = (-150 : -50); %for "focused" plot

sdfAcc = NaN(NUM_CELLS, N_SAMP);
sdfFast = NaN(NUM_CELLS, N_SAMP);

for cc = 1:NUM_CELLS
  fprintf('%s - %s\n', ninfo(cc).sess, ninfo(cc).unit)
  
  kk = ismember({binfo.session}, ninfo(cc).sess);
  
  %compute single-trial SDF
  SDFcc = compute_spike_density_fxn(spikes(cc).SAT);
  
  %index by isolation quality
  idxIso = identify_trials_poor_isolation_SAT(ninfo(cc), binfo(kk).num_trials, 'task','SAT');
  %index by trial outcome
%   idxCorr = ~(binfo(kk).err_dir | binfo(kk).err_time | binfo(kk).err_hold | binfo(kk).err_nosacc);
  idxCorr = ~(binfo(kk).err_time | binfo(kk).err_hold | binfo(kk).err_nosacc);
  %index by condition
  idxAcc = ((binfo(kk).condition == 1) & idxCorr & ~idxIso);
  idxFast = ((binfo(kk).condition == 3) & idxCorr & ~idxIso);
  %index by response direction relative to RF and MF
  [visField,~] = determineFieldsVisMove( ninfo(cc) );
  idxRF = ismember(moves(kk).octant, visField);
  
  %split single-trial SDF by condition
  sdfAccST = SDFcc(idxAcc & idxRF, T_STIM);
  sdfFastST = SDFcc(idxFast & idxRF, T_STIM);
  
  %compute mean SDF
  sdfAcc(cc,:) = mean(sdfAccST);
  sdfFast(cc,:) = mean(sdfFastST);
  
end%for:cells(cc)

%% Plotting

%normalization
sdfAcc = sdfAcc ./ [nstats.NormFactor_All]';
sdfFast = sdfFast ./ [nstats.NormFactor_All]';

%split neurons by level of search efficiency
ccMore = ([ninfo.taskType] == 1);   NUM_MORE = sum(ccMore);
ccLess = ([ninfo.taskType] == 2);   NUM_LESS = sum(ccLess);
sdfAccMore = sdfAcc(ccMore,:);     sdfFastMore = sdfFast(ccMore,:);
sdfAccLess = sdfAcc(ccLess,:);     sdfFastLess = sdfFast(ccLess,:);

T_STIM = T_STIM - 3500;

%time from stimulus for plotting close-up of baseline
IDX_FOCUSED = ismember(T_STIM, T_FOCUSED);

%compute common y-axis scale
tmp = [mean(sdfAccMore) mean(sdfFastMore) mean(sdfAccLess) mean(sdfFastLess)];
yLim = [(min(tmp)-0.05) , (max(tmp)+0.05)];

figure()

%Less efficient
subplot(2,2,1); hold on %from stimulus
plot([0 0], yLim, 'k:')
shaded_error_bar(T_STIM, mean(sdfAccLess), std(sdfAccLess)/sqrt(NUM_LESS), {'r-', 'LineWidth',1.25})
shaded_error_bar(T_STIM, nanmean(sdfFastLess), nanstd(sdfFastLess)/sqrt(NUM_LESS), {'-', 'Color',[0 .7 0], 'LineWidth',1.25})
xlabel('Time from array (ms)'); ytickformat('%2.1f')

subplot(2,2,2); hold on %focused look at baseline
plot(T_FOCUSED, mean(sdfAccLess(:,IDX_FOCUSED)), 'r-', 'LineWidth',1.25)
plot(T_FOCUSED, nanmean(sdfFastLess(:,IDX_FOCUSED)), '-', 'Color',[0 .7 0], 'LineWidth',1.25)

%More efficient
subplot(2,2,3); hold on %from stimulus
plot([0 0], yLim, 'k:')
shaded_error_bar(T_STIM, mean(sdfAccMore), std(sdfAccMore)/sqrt(NUM_MORE), {'r-', 'LineWidth',0.75})
shaded_error_bar(T_STIM, mean(sdfFastMore), std(sdfFastMore)/sqrt(NUM_MORE), {'-', 'Color',[0 .7 0], 'LineWidth',0.75})
ylabel('Norm. activity'); ytickformat('%2.1f')

subplot(2,2,4); hold on %focused look at baseline
plot(T_FOCUSED, mean(sdfAccMore(:,IDX_FOCUSED)), 'r-', 'LineWidth',0.75)
plot(T_FOCUSED, mean(sdfFastMore(:,IDX_FOCUSED)), '-', 'Color',[0 .7 0], 'LineWidth',0.75)

ppretty([8,1.8])

end%fxn:plotBaselineSDF_SAT()


function [visField , moveField] = determineFieldsVisMove( ninfo )

if (isempty(ninfo.visField) || ismember(9, ninfo.visField)) %non-specific RF
  visField = (1:8);
else %specific RF
  visField = ninfo.visField;
end

if (isempty(ninfo.moveField) || ismember(9, ninfo.moveField)) %non-specific MF
  moveField = (1:8);
else %specific MF
  moveField = ninfo.moveField;
end

end%util:determineFieldsVisMove()
