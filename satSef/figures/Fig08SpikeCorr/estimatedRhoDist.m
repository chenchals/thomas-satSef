% plot observed rho on the estimted rho distribution
% estimated rohs are for 1000 sub sampling of chosing
% 1. 40 random trials or
% 2. 80 random trials
% from set fo trials for a given sat-outcome combination.

% variables from spk corr are:
% Pair_UID
% X_monkey
% X_area
% Y_area
% condition
% nTrials
% rhoRaw : observed Rho
% subSamplIdxs_nTrials_40 : 1000 sets of 40 random-choices 
% rhoVecSubSampl_nTrials_40 : vector of subSampled rho (40 trials)
% rhoEstRaw_nTrials_40 : mean rho estimated from rhoVecSubSampl_nTrials_40 
% ci95_nTrials_40 : 95% CI from rhoVecSubSampl_nTrials_40 for rhoRaw             
% rhoRawInCi95_nTrials_40 : observed rhoRaw within CI ci95_nTrials_40
% subSamplIdxs_nTrials_80 : 1000 sets of 80 random-choices 
% rhoVecSubSampl_nTrials_80 :  vector of subSampled rho (80 trials)
% rhoEstRaw_nTrials_80 : mean rho estimated from rhoVecSubSampl_nTrials_80 
% ci95_nTrials_80 : 95% CI from rhoVecSubSampl_nTrials_80 for rhoRaw             
% rhoRawInCi95_nTrials_80 : observed rhoRaw within CI ci95_nTrials_80
%%
spkCorr = load('dataProcessed/satSefPaper/rscSubSampl1K_PostSaccade.mat');
spkCorr = spkCorr.spkCorr;
spkCorr = spkCorr(:,{...
    'Pair_UID'
    'X_monkey'
    'X_area'
    'Y_area'
    'condition'
    'nTrials'
    'rhoRaw'
    'signifRaw_05'
    'subSamplIdxs_nTrials_40'
    'rhoVecSubSampl_nTrials_40'
    'rhoEstRaw_nTrials_40'
    'rhoEstSem_nTrials_40'
    'ci95_nTrials_40'
    'rhoRawInCi95_nTrials_40'
    'subSamplIdxs_nTrials_80'
    'rhoVecSubSampl_nTrials_80'
    'rhoEstRaw_nTrials_80'
    'rhoEstSem_nTrials_80'
    'ci95_nTrials_80'
    'rhoRawInCi95_nTrials_80'});

%%
oPdfDir = 'estimatedRhoDistrib';
if ~exist(oPdfDir,'dir')
    mkdir(oPdfDir)
end

satConds = {'Accurate','Fast'};
outcomes = {'Correct','ErrorChoice','ErrorTiming'};
crossPairs = unique(spkCorr.Pair_UID);
for cp = 1:numel(crossPairs) %10
    pairUid = crossPairs{cp};
    oPdfFile = fullfile(oPdfDir,[pairUid '_estimatedRscDistrib.pdf']);
    H_plots = figTemplate();
    pl = 0;
    for o = 1:numel(outcomes)
       outcome = outcomes{o};
       pl = o;
       labelOutcome = true;
       if o == 1
           labelSat = true;
       end
       for s = 1: numel(satConds)
        satCond = satConds{s};
            satOutcome = [satCond outcome];
            idx = ismember(spkCorr.Pair_UID,pairUid) & ismember(spkCorr.condition,satOutcome);
            
            obs = spkCorr.rhoRaw(idx);
            est40 = spkCorr.rhoEstRaw_nTrials_40(idx);
            est80 = spkCorr.rhoEstRaw_nTrials_80(idx);
            ci40 = spkCorr.ci95_nTrials_40{idx};
            ci80 = spkCorr.ci95_nTrials_80{idx};
            
            trlIdxEst40 = spkCorr.subSamplIdxs_nTrials_40{idx};
            trlIdxEst80 = spkCorr.subSamplIdxs_nTrials_80{idx};
            rawEst40 = spkCorr.rhoVecSubSampl_nTrials_40{idx};
            rawEst80 = spkCorr.rhoVecSubSampl_nTrials_80{idx};
            
            sig = spkCorr.signifRaw_05(idx);
            % histogram for random selection of 40/80 trials
            trialBins = 1:spkCorr.nTrials(idx);
            rhoBins = -0.5:0.01:0.5;
            cTrials40 = histcounts(trlIdxEst40,'BinEdges',[0.5 trialBins+0.5]);
            cRhoEst40 = histcounts(rawEst40,'BinEdges',[rhoBins(1)-0.005 rhoBins+0.005]);
            cTrials80 = histcounts(trlIdxEst80,'BinEdges',[0.5 trialBins+0.5]);
            cRhoEst80 = histcounts(rawEst80,'BinEdges',[rhoBins(1)-0.005 rhoBins+0.005]);
            
            % plot distribution of randomly selected trial index
            axes(H_plots(pl))
            bar(trialBins,[cTrials40;cTrials80]','BarWidth',1)
            xlabel('Trial Idx');
            ylabel('Count of trial Idx')
           if labelOutcome
                title(outcome,'FontWeight','bold');
                labelOutcome = false;
            end
            if labelSat
                pos = get(gca,'Position');
                axes('Position',[0.001 pos(2)-0.01 0 0])
                text(0,0,upper(satCond),'Rotation',90,'FontWeight','bold',...
                'FontAngle','italic','FontSize',12,...
                'HorizontalAlignment','center','VerticalAlignment','top')
            end
            % plot estimated rho distribution
            % offset by 3 to next plot in the column
            pl = pl + 3;
            axes(H_plots(pl))
            bar(rhoBins,cRhoEst40,'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off')
            hold on
            bar(rhoBins,cRhoEst80,'FaceAlpha',0.4,'EdgeColor','none','HandleVisibility','off')
            xlabel('Estimated Rsc');
            ylabel('count of Rsc')
            
            % Annotate rho-distribution
            yLim = get(gca,'YLim');
            yMin = yLim(1);
            yMax = yLim(2);
            yLim110 = [yMin yMax*1.1];
            
            legTxt = {};
            
            plot(obs,yMax*0.9,'*r','MarkerSize',10)
            legTxt = [legTxt sprintf('\\mu_{obs.} %+0.3f',obs)];
            plot([est40 est40],yLim,'-b','LineWidth',0.75)
            legTxt = [legTxt sprintf('\\mu_{est40.} %+0.3f',est40)];
            plot([est80 est80],yLim,'--r','LineWidth',1)
            legTxt = [legTxt sprintf('\\mu_{est80.} %+0.3f',est80)];
            plot(ci40,[yMax yMax].*1.08,'-b','LineWidth',1.5)
            legTxt = [legTxt ['CI_{est40.} ' sprintf('%+0.3f ',ci40)]];
            plot(ci80,[yMax yMax].*1.04,'--r','LineWidth',2)
            legTxt = [legTxt ['CI_{est40.} ' sprintf('%+0.3f ',ci80)]];

            set(gca,'YLim',yLim110);             
            [~,objH] = legend(legTxt,'Location','northeast','Interpreter','tex');
            objhl = findobj(objH, 'type', 'line'); %// objects of legend of type line
            set(objhl, 'Markersize', 8); %// set marker size as desired
            
            hold off
            % offset by 3 to next plot in the column 
            pl = pl + 3;
        end
    end
   ha = annotation('textbox','String',pairUid,...
        'Position',[0.02 0.97 0.01 0.01],'LineStyle','none',...
        'FontSize',12,'FontWeight','bold','Interpreter','none');
    set(gcf,'Position',[120 120 1400 900]);
    set(gcf,'PaperOrientation','landscape')
    drawnow
    saveFigPdf(oPdfFile);
    delete(gcf);
end


%%
function [hA] = figTemplate()
%%
figure;
hA = tight_subplot(4, 3, [.05 .05],[0.05 0.08],[0.05 0.05]);
for ii = 1:numel(hA)
    axes(hA(ii))
    title(sprintf('Plot %d',ii));
end
end




