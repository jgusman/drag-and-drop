%% Gesture Hero LDA Over Time
% generates main figures 2(B) and 2(C) and supplemental figures S2(B), S2(C), S3(B), S3(C), S4(B), and S4(C) 

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------

% set figure save subfolder
ddFigHelper.ResetParams()
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'LDAOverTime'))

%% Set Common Sweep Params
sweepInfo.attemptWindow   = [-24 0]; % What goes into the lda for decoding. [-24 0] is 0.5 seconds back.
sweepInfo.attemptOffset   = [-50 85]; % Sweep boundaries relative to [attempt-trial-start, attempt-trial-stop].
sweepInfo.noAttemptWindow = [25 150]; % No-attempt period relative to noActionStartStops start. (attempt trial end).
sweepInfo.stepSize = 5; % <--- 5 (number of 20ms time steps to move in each step of calculation, increase for faster, crude calculation)
sweepInfo.featInds = 1:size(sesData_GH_T11(1).feat,2);%1:384;
sweepInfo.enforceTransitionMatrix = 0; % 0 = No HMM. [] = Yes HMM. Transition matrix from train data trial transitions.


allFeatNames = {'ncTX','spikePower','STFT'};
featTypeInds = {1:192,193:384,385:1344};



%% FIGURE 2(B) and 2(C)
% T11 - Calculate LDA Over Time Gesture and Attempt - TX and SP

featNames = {'ncTX','spikePower'};
selFeat = [featTypeInds{ismember(allFeatNames,featNames)}];
sweepPerSesT11 = ComparePerformanceOverTime_SelFeat(sesData_GH_T11,sweepInfo,selFeat);

% Plot gesture and attempt separately
[fh1, fh2] = PlotGestureAndAttemptOverTimeSeparatePlots(sweepPerSesT11);

% Save figures
ddFigHelper.SaveFigure(fh1,sprintf('Figure 2(B) - Sliding LDA - Gesture with CI (T11), %s',[featNames{:}]))
ddFigHelper.SaveFigure(fh2,sprintf('Figure 2(C) - Sliding LDA - Attempt with CI (T11), %s',[featNames{:}]))



%% FIGURE S2(B) and S2(C)
% T11 - Calculate LDA Over Time Gesture and Attempt *BY FEATURE TYPE*

ddFigHelper.LogPrints('PerfOverTimePeaks_EachFeature_T11_7g.txt')
sweepPerSesPerFeatTypeT11 = ComparePerformanceOverTimeByFeat(sesData_GH_T11,sweepInfo);
ddFigHelper.LogPrints('off')

% Plot by Feature Type - 7 Gestures - 4s duration trials
PlotLDAOverTimeByFeature(sweepPerSesPerFeatTypeT11,'gest');  
legend off
ddFigHelper.SaveFigure(sprintf('Figure S2(B) - Sliding LDA - Gesture (T11), SeparateFeats (%s) 7g', [allFeatNames{:}]))

PlotLDAOverTimeByFeature(sweepPerSesPerFeatTypeT11,'attempt');
legend off
ddFigHelper.SaveFigure(sprintf('Figure S2(C) - Sliding LDA - Attempt (T11), SeparateFeats (%s) 7g', [allFeatNames{:}]))



%% FIGURE S3(B) and S3(C)
% T5 - Calculate LDA Over Time Gesture and Attempt - TX and SP

featNames = {'ncTX','spikePower'};
selFeat = [featTypeInds{ismember(allFeatNames,featNames)}];
sweepPerSesT5 = ComparePerformanceOverTime_SelFeat(sesData_GH_T5,sweepInfo,selFeat);

%% Plot gesture and attempt separately
[fh1, fh2] = PlotGestureAndAttemptOverTimeSeparatePlots(sweepPerSesT5);

% Save figures
ddFigHelper.SaveFigure(fh1,sprintf('Figure S3(B) - Sliding LDA - Gesture (T5), %s',[featNames{:}]))
ddFigHelper.SaveFigure(fh2,sprintf('Figure S3(C) - Sliding LDA - Attempt (T5), %s',[featNames{:}]))



%% FIGURE S4(B) and S4(C) 
% T11 - Calculate LDA Over Time Gesture and Attempt - TX and SP - 3 GESTURES

gestures = {'Right index finger - down', 'Right hand - ok', 'Right hand - power grasp'};

featNames = {'ncTX','spikePower'};
selFeat = [featTypeInds{ismember(allFeatNames,featNames)}];
sweepPerSesT11_3g = ComparePerformanceOverTime_SelFeat(sesData_GH_T11,sweepInfo,selFeat,gestures);

% Plot gesture and attempt separately
[fh1, fh2] = PlotGestureAndAttemptOverTimeSeparatePlots(sweepPerSesT11_3g);

% Save figures
ddFigHelper.SaveFigure(fh1,sprintf('Figure S4(B) - Sliding LDA - Gesture (T11 - 3 gest), %s',[featNames{:}]))
ddFigHelper.SaveFigure(fh2,sprintf('Figure S4(C) - Sliding LDA - Attempt (T11 - 3 gest), %s',[featNames{:}]))





%% %%%% HELPER FUNCTIONS %%%% %%

function sweepPerSes = ComparePerformanceOverTime_SelFeat(sesData,sweepInfo,selFeat,gestures)

%% Gesture and Attempt performance over time.
% Sweep LDA performance across attempt period of the trial.
% For each attempt duration,
% Sweep for gesture labels
% Sweep for Atempt/no attempt labels

sweepInfo.featInds = selFeat;

if nargin < 4
    uGestLabels = unique(sesData.taskInfo.labels);
    gestures = uGestLabels;
end

infoStr = 'No HMM. No attempt is 0.5-3 sec after attempt offset';

tt = tic;
for sesI = 1:length(sesData)
    fprintf('\n\nSweeping %s Session %d\n', upper(sesData(sesI).participant), sesData(sesI).sessionNum);
    
    sweepInfo.labelType = 'gestures';
    sweepGest = ComparePerformanceOverTime(sesData(sesI),sweepInfo,gestures);

    sweepInfo.labelType = 'attempt';
    sweepAttempt = ComparePerformanceOverTime(sesData(sesI),sweepInfo,gestures);
    
    sweepPerSes(sesI).infoStr = infoStr;
    sweepPerSes(sesI).gesturesUsed = gestures;
    sweepPerSes(sesI).gest = sweepGest;
    sweepPerSes(sesI).attempt = sweepAttempt;
    sweepPerSes(sesI).sweepInfo = sweepInfo;
end
fprintf('\nDone!\nSweep took %.1f seconds (%.1f min.)\n', toc(tt),toc(tt)/60);

end


function [fh1,fh2] = PlotGestureAndAttemptOverTimeSeparatePlots(plotSwp)

    % Setup figure
    fwh = [3.2 1.7];

    % Axes options
    axArgs = ddFigHelper.GetDefaultAxArgs();
    axArgs.marg_w = [0.13 0.05]; % Update axes width boundaries
    axArgs.marg_h = [0.2 0.08]; % Update axes height boundaries
    axArgs.xTickOffset = -1; % Make x ticks closer to axes, I think.
    
    % Plot options
    plotOpt.mainM = '-';
    plotOpt.showCI = '--';
    plotOpt.showLegend = 1;
    

    fh1 = ddFigHelper.CreateFigure(32240,fwh);
    clf
    axs1 = ttight_subplot(1,1,axArgs);
    
    if length(plotSwp(1).gest.perf{1}(1).accuracy)==8 % 7 gesture plot for fig 1/2
        plotOpt.forceYLim = [10 80];
    else  %i.e. if 3 gestures
        plotOpt.forceYLim = [20 100];
    end

    plotOpt.clrs = ddHelper.durationColors;
    
    % Plot Gesture sweep
    fprintf('Gesture Over Time\n')
    bb = [plotSwp.gest];
    PlotSweepOverTimeMultiDuration(bb,axs1(1,:),plotOpt);


    % Plot Attempt sweep
    fh2 = ddFigHelper.CreateFigure(32241,fwh);
    clf
    axs2 = ttight_subplot(1,1,axArgs);

    plotOpt.forceYLim = [40 100];
    
    plotOpt.clrs = ddHelper.attemptDurationColors;

    fprintf('Attempt Over Time\n')
    bb = [plotSwp.attempt];
    PlotSweepOverTimeMultiDuration(bb,axs2(1,:),plotOpt);


    % Hacky clean up
    axs1.XLim = [-1 6];
    axs1.XTick = RowColon(axs1.XLim);
    axs1.XLabel.Units = 'normalized'; %idk if necessary?
    axs1.XLabel.FontSize = 10;
    axs2.XLim = [-1 6];
    axs2.XTick = RowColon(axs2.XLim);
    axs2.XLabel.Units = 'normalized';
    axs2.XLabel.FontSize = 10;

end


function sweepPerSesPerFeatType = ComparePerformanceOverTimeByFeat(sesData,sweepInfo_FT,gestures)
    
    if nargin < 3
        uGestLabels = unique(sesData(1).taskInfo.labels);
        gestures = uGestLabels;
    end

    numFeatTypes = length(sweepInfo_FT.featInds)/192;
    featTypeNames = {'TX','SP','0-11Hz','12-19Hz','20-39Hz','40-128Hz','129-250Hz'};
                        
    tt = tic;
        
    for featTypeI = 1:numFeatTypes
        for sesI = 1:length(sesData)
        
            fprintf('\n\nSweeping %s Session %d\n', upper(sesData(sesI).participant), sesData(sesI).sessionNum);

            sweepInfo_FT.featInds = ((featTypeI-1)*192+1) : (featTypeI*192);
    
            sweepInfo_FT.labelType = 'gestures';
            sweepGest = ComparePerformanceOverTime(sesData(sesI),sweepInfo_FT,gestures);
        
            sweepInfo_FT.labelType = 'attempt';
            sweepAttempt = ComparePerformanceOverTime(sesData(sesI),sweepInfo_FT,gestures);
            
            infoStr = 'No HMM. No attempt is 0.5-3 sec after attempt offset.';
        
            sweepPerSesPerFeatType(featTypeI,sesI).infoStr = infoStr;
            sweepPerSesPerFeatType(featTypeI,sesI).featType = featTypeNames{featTypeI};
            sweepPerSesPerFeatType(featTypeI,sesI).gesturesUsed = gestures;
            sweepPerSesPerFeatType(featTypeI,sesI).gest = sweepGest;
            sweepPerSesPerFeatType(featTypeI,sesI).attempt = sweepAttempt;
            sweepPerSesPerFeatType(featTypeI,sesI).sweepInfo = sweepInfo_FT;
        end
    
        [fh1, fh2] = PlotGestureAndAttemptOverTimeSeparatePlots(sweepPerSesPerFeatType(featTypeI,:));
        
        legend(fh1.Children(2),'off') %remove legend
        legend(fh2.Children(2),'off') %remove legend
    
        % Save figures
%         ddFigHelper.SaveFigure(fh1,sprintf('Sliding LDA - Gesture with CI (%s), %s, %dg',upper(sesData(1).participant),featTypeNames{featTypeI},length(gestures)))
%         ddFigHelper.SaveFigure(fh2,sprintf('Sliding LDA - Attempt with CI (%s), %s, %dg',upper(sesData(1).participant),featTypeNames{featTypeI},length(gestures)))
%     
    end
    fprintf('\nDone!\nSweep took %.1f seconds (%.1f min.)\n', toc(tt),toc(tt)/60);
        
end


function fh = PlotLDAOverTimeByFeature(plotSwp,gestORattempt)

    % Setup figure
    fwh = [3.2 1.7];

    % Axes options
    axArgs = ddFigHelper.GetDefaultAxArgs();
    axArgs.marg_w = [0.13 0.05]; % Update axes width boundaries
    axArgs.marg_h = [0.2 0.08]; % Update axes height boundaries
    axArgs.xTickOffset = -1; % Make x ticks closer to axes, I think.
    
    % Plot options
    plotOpt.mainM = '-';
    plotOpt.showCI = '--';
    plotOpt.showLegend = 1;
    plotOpt.lw = 1;

    try
        fig = gcf;
        figNum = fig.Number+1;
    catch
        figNum = 78901;
    end

    fh = ddFigHelper.CreateFigure(figNum,fwh);
    clf
    axs = ttight_subplot(1,1,axArgs);
    
    
    if strcmp(gestORattempt,'attempt')
        plotOpt.forceYLim = [40 100];
        stopLineColor = ddHelper.attemptDurationColors(3,:);
    else
        if length(plotSwp(1).gest.perf{1}(1).accuracy)==8 % 7 gesture plot for fig 1/2
            plotOpt.forceYLim = [10 80];
        else  %i.e. if 3 gestures
            plotOpt.forceYLim = [20 100];
        end
        stopLineColor = ddHelper.durationColors(3,:);
    end
    
    clrVar = ['featTypeColors_',gestORattempt];
    clrVar = 'featTypeColors'; %same
    plotOpt.clrs = ddHelper.(clrVar);
    
    selDur = 200; 
    for i = 1:size(plotSwp,1)
        for j = 1:size(plotSwp,2)
            swp = [plotSwp(i,j).(gestORattempt)];
            durInd = find(ismember(swp.durations,selDur));
            swp.durations = swp.durations(durInd);
            swp.offsets = swp.offsets(durInd);
            swp.perf = swp.perf(durInd);
            swp_all(i,j) = swp;
        end
    end
    
    featPltOrder = [1 2 7 6 5 4 3];
    swp_all = swp_all(featPltOrder,:);
    
    ax = PlotSweepOverTimeMultiDuration(swp_all,axs(1,:),plotOpt);

    ax.Children(13).Color = stopLineColor;
    ax.Children(14).MarkerFaceColor = stopLineColor;
    
    lgStr = {'','','',ddHelper.featTypeNames{featPltOrder}};
    legend(lgStr,'FontSize',5,'Location','best')


    % Hacky clean up
    axs.XLim = [-1 6];
    axs.XTick = RowColon(axs.XLim);
    axs.XLabel.Units = 'normalized'; %idk if necessary?
    axs.XLabel.FontSize = 10;


end




