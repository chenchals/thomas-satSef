function [] = figure1()
%FIGURE1 Plot Fig1b,Fig1c
% needs Behavior data
% see also FIG01B_PLOT_ERRRATE_X_RT, 
%          FIG01C_PLOT_BEHAV_X_TRIAL, 
%          GETOPT,
%          PPRETTY
% Needs behavior data
%%
binfoFile = 'dataProcessed/dataset/dataBehavior_SAT.mat';
binfoAll = load(binfoFile);
binfo = binfoAll.binfoSAT;
pSacc = binfoAll.primarySaccade;
sSacc = binfoAll.secondSaccade;

%% Call Error Rate by Reaction time plot
Fig01B_Plot_ErrRate_X_RT(binfo,pSacc);

%%
Fig01C_Plot_Behav_X_Trial(binfo,pSacc);


end

