%% Gesture and Attempt Selectivity Over Time
% generates Figs 2A, 3A-B, S2A, S3A, S4A, S5

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------

% set figure save subfolder
ddFigHelper.ResetParams()
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'SelectivityOverTime'))

%% Set Common Vars

winOnsetOffset = [-50 85]; % [-200 100] %(steps) test onset, test offset from each hold duration
testWin = -14:0; % (time steps) averaged over before testing 
stepSize = 1; % <--- 1 (number of 20ms time steps to move in each step of calculation, increase for faster, crude calculation)
sigAlpha = 0.001; % Sig threshold
figWH = [3.2 1.6]; % Default fig [width height] in inches


%% Calculate Selective Features for T11 (all 7 gestures and 3 gesture subset) and T5 (all 3 gestures)

for ii = 1:length(sesData_GH_T11) % 2 sessions

% T11 Selectivity Over Time - all 7 gestures
gestures = unique(sesData_GH_T11(1).taskInfo.labels);
gestSel_T11_7g{ii} = CalcSelectiveOverTime(sesData_GH_T11(ii),winOnsetOffset,testWin,stepSize,false,false,gestures);   % gesture selectivity
attemptSel_T11_7g{ii} = CalcSelectiveOverTime(sesData_GH_T11(ii),winOnsetOffset,testWin,stepSize,true,true,gestures);  % attempt selectivity

% T11 Selectivity Over Time - just 3 gestures
gestures =  {'Right index finger - down', 'Right hand - ok', 'Right hand - power grasp'};
gestSel_T11_3g{ii} = CalcSelectiveOverTime(sesData_GH_T11(ii),winOnsetOffset,testWin,stepSize,false,false,gestures);    % gesture selectivity
attemptSel_T11_3g{ii} = CalcSelectiveOverTime(sesData_GH_T11(ii),winOnsetOffset,testWin,stepSize,true,true,gestures);   % attempt selectivity

end

for ii = 1:length(sesData_GH_T5)

% T5 Selectivity Over Time - all 3 gestures
gestures = unique(sesData_GH_T5(1).taskInfo.labels);
gestSel_T5_3g{1} = CalcSelectiveOverTime(sesData_GH_T5(ii),winOnsetOffset,testWin,stepSize,false,false,gestures);   % gesture selectivity
attemptSel_T5_3g{1} = CalcSelectiveOverTime(sesData_GH_T5(ii),winOnsetOffset,testWin,stepSize,true,true,gestures);  % attempt selectivity

end

%% Plot Main Paper Figures

% Figure 2(A) - T11 Gesture Selectivity Over Time - 7 gestures - 1s 2s 4s - just TX and SP feats
PlotKW_selFeats(gestSel_T11_7g,sigAlpha,[],figWH,1:384);
set(gcf, 'Name', 'Figure 2(A)', 'NumberTitle', 'off');
ddFigHelper.SaveFigure('Figure 2(A) - T11 Gesture Selectivity Over Time - 7 gestures - 1s 2s 4s - just TX and SP feats')

% Figure 3(A) - T11 Attempt Selectivity Over Time - 3 gestures - 1s 2s 4s - just TX and SP feats
PlotKW_selFeats(attemptSel_T11_3g,sigAlpha,[],figWH,1:384);
set(gcf, 'Name', 'Figure 3(A)', 'NumberTitle', 'off');

% Figure 3(B) - T5 Attempt Selectivity Over Time - 3 gestures - 1s 2s 4s - just TX and SP feats
PlotKW_selFeats(attemptSel_T5_3g,sigAlpha,[],figWH,1:384);
set(gcf, 'Name', 'Figure 3(B)', 'NumberTitle', 'off');


%% Plot Supplemental Figures 

% Figure S2(A) - T11 Gestures Selectivity Over Time - 7 gestures - 4s only - by feature type
PlotKW_ByFeatType_SamePlot(gestSel_T11_7g,sigAlpha,figWH)
set(gcf, 'Name', 'Figure S2(A)', 'NumberTitle', 'off');

% Figure S3(A) - T11 Gesture Selectivity Over Time - 3 gestures - 1s 2s 4s - just TX and SP feats
PlotKW_selFeats(gestSel_T11_3g,sigAlpha,[],figWH,1:384);
set(gcf, 'Name', 'Figure S3(A)', 'NumberTitle', 'off');

% Figure S4(A) - T5 Gesture Selectivity Over Time - 3 gestures - 1s 2s 4s - just TX and SP feats
PlotKW_selFeats(gestSel_T5_3g,sigAlpha,[],figWH,1:384);
set(gcf, 'Name', 'Figure S4(A)', 'NumberTitle', 'off');

% Figure S5(A) - T11 Attempt Selectivity Over Time - 7 gestures - 4s only - by feature type
PlotKW_ByFeatType_SamePlot(attemptSel_T11_7g,sigAlpha,figWH)
set(gcf, 'Name', 'Figure S5(A)', 'NumberTitle', 'off');

% Figure S5(B) - T11 Attempt Selectivity Over Time - 3 gestures - 4s only - by feature type
PlotKW_ByFeatType_SamePlot(attemptSel_T11_3g,sigAlpha,figWH)
set(gcf, 'Name', 'Figure S5(B)', 'NumberTitle', 'off');

% Figure S5(C) - T5 Attempt Selectivity Over Time - 3 gestures - 4s only - by feature type
PlotKW_ByFeatType_SamePlot(attemptSel_T5_3g,sigAlpha,figWH)
set(gcf, 'Name', 'Figure S5(C)', 'NumberTitle', 'off');









%% %%%% HELPER FUNCTIONS %%%% %%


function kwPerDuration = CalcSelectiveOverTime(sesData, winOnsetOffset, testWin, stepSize,includeBL,labelAllAttempt,gestures)
% kwPerDuration = CalcSelective(sesData, winOnsetOffset, stepSize)
% 
% Compare sliding Kruskal Wallis (KW) test across all neural features for each
% unique attempt duration. 
% 
% Inputs
%   sesData
%       Struct with the following substructs
%       .feat - neural features
%       .taskInfo - taskInfo struct with 
%           .startStops
%           .labels
%           .trialDurationRounded
%           .prctNS5Outliers
%   
%   winOnsetOffset length: 2
%       relative [onset, offset] steps around events to sweep from attempt
%       start to attempt stop.
% 
%       example: [-50 50] will sweep from 50 time steps before the trial
%       starts to 50 time steps after the respective trial attempt
%       duration.
% 
%   testWin (default: -14:0, which is 300ms for 20ms steps)
%       Window averaged over when computing the statistical test.
% 
%   stepSize (default: 1)
%       How many steps to take between each sliding onset/offset window.
%          
%   gestures (default: unique(taskInfo.labels)
%       What gestures to include in selectivity calculation
%
% Output
%   kwPerDuration
%       cell array of statistical test structures. One entry for each
%       unique duration.
%       .p contains the p values for the statistcal test.

if nargin < 3
    testWin = -14:0;
end
if nargin < 4
    stepSize = 1;
end
if nargin < 5
    includeBL = false;
end
if nargin < 6
    labelAllAttempt = false;
elseif isempty(labelAllAttempt)
    labelAllAttempt = false;
end
if nargin < 7
    uGestLabels = unique(sesData.taskInfo.labels);
    gestures = uGestLabels;
end


outlierThreshold = 3;

taskInfo = sesData.taskInfo;
durs = taskInfo.trialDurationRounded;
[cnt_uDurations,uDurations] = hist(durs,unique(durs));

uDurations(cnt_uDurations<3) = [];  %remove one-off durations (or at least those with less than 3 trials (arbitrary) )



featInds = 1:size(sesData.feat,2);
for jj = 1:length(uDurations)

    analyzeWindow = winOnsetOffset(1):stepSize:(uDurations(jj)+winOnsetOffset(2));

    % Get average outlier per trial
    % Exclude
    outlierWin = analyzeWindow(1):analyzeWindow(end);
    outlierPerTrial = mean(GetEpochOfData2( taskInfo.prctNS5Outliers, taskInfo.startStops(:,1), outlierWin,'edt' ),3);
    excludeTrials = outlierPerTrial(:) > outlierThreshold;
    
    % Select all trials for a given duration
    selDur = ismember(taskInfo.trialDurationRounded, uDurations(jj));

    % Select trials only certain gestures
    selGest = ismember(taskInfo.labels,gestures);
    
    selTrl = selDur & selGest & ~excludeTrials;
    
    fprintf('Percent of Trials Excluded: %0.0f%%\n', mean(excludeTrials)*100)
    
    % Set up Kruskal Wallis stat test for each feature.
    %
    % No baseline events means we compare only attempt periods, i.e. a
    % selectivity (not responsiveness) test.
    %
    % Test each neural feature's average response across relative
    % window, winInds. Repeat for each slide window value in
    % analyzeWindow.

    if includeBL 
        kwOpt.baselineStatWinInds = -14:0;
        % Get average outlier per NA trial
        outlierWin_BL = kwOpt.baselineStatWinInds;
        outlierPerTrial_BL = mean(GetEpochOfData2( taskInfo.prctNS5Outliers, taskInfo.noActionStartStops(:,1), outlierWin_BL,'edt' ),3);
        excludeTrials_BL = outlierPerTrial_BL(:) > outlierThreshold;
        selTrl_BL = ~selDur & selGest & ~excludeTrials_BL; %using no action epochs from other durations to avoid overlap
        kwOpt.baselineEvents = taskInfo.noActionStartStops(selTrl_BL,1) + 65; % using default baselineStatWinInds (-14:0), uses window 1.0-1.3s after attempt ends
    end
    
    kwOpt.events = taskInfo.startStops(selTrl,1);
    kwOpt.labels = removecats(taskInfo.labels(selTrl));
    kwOpt.winInds = testWin; % What to average over
    kwOpt.slideWinInds = analyzeWindow; % Relative points to analyze around events
    kwOpt.plotOn = false;
    
    if labelAllAttempt
        kwOpt.labels(:) = categorical({'Attempt'});
    end

    [kwP,kwOpt] = ComputeSignificantChAndPlotPeriStim(sesData.feat(:,featInds), kwOpt);
    
    % kwP is size: [nSlideWin x nFeatures]
    kwOpt.p = kwP;
    
    kwPerDuration{jj} = kwOpt;
end

end

function ax = PlotKW_selFeats(kwData,sigAlpha,figNum, figWH,selFeats)

    plotXLim = [-1 6]; % seconds rel to go cue

    kws = {};
    for s = 1:length(kwData)
        data = kwData{s};
        for d = 1:length(data)
            data{d}.p = data{d}.p(:,selFeats); %reorganize feat order for plotting purposes
        end
        kws{s} = data;
    end

    if length(unique(kws{1}{1}.labels)) == 1
        durClrs = ddHelper.attemptDurationColors;
    else
        durClrs = ddHelper.durationColors;
    end
    
    ax = PlotKW(kws,sigAlpha,figNum, figWH, durClrs);

    ax.XLim = plotXLim; 

end


function PlotKW_ByFeatType_SamePlot(kwPerDurationPerSes,sigAlpha,figWH)

    figNum = [];
    
    selDur = [200]; %4sec trials
    durInd = find(ismember([50 100 200],selDur));
    featPltOrder = [1 2 7 6 5 4 3];
    
    newFeatInds = [];
    for f = featPltOrder
        featInds = ((f-1)*192+1) : (f*192);
        newFeatInds = [newFeatInds,featInds];
    end
    
    kws = {};
    for s = 1:length(kwPerDurationPerSes)
        data = kwPerDurationPerSes{s}(durInd);
        for d = 1:length(durInd)
            data{d}.p = data{d}.p(:,newFeatInds); %reorganize feat order for plotting purposes
        end
        kws{s} = data;
    end
    
    ax = PlotKW(kws,sigAlpha,figNum, figWH, ddHelper.featTypeColors,true,selDur);
    
    if length(unique(kws{1}{1}.labels)) == 1
        stopLineColor = ddHelper.attemptDurationColors(3,:);
    else
        stopLineColor = ddHelper.durationColors(3,:);
    end

    ax.Children(9).Color = stopLineColor;
    ax.Children(10).MarkerFaceColor = stopLineColor;
    
    xlim([-1 6])
    
    lgStr = {'','','',ddHelper.featTypeNames{featPltOrder}};
    legend(lgStr,'FontSize',5,'Location','best')

end

