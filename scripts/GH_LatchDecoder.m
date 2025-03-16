%% Latch Decoder example trial and comparison with Gesture Decoder (Figs 4B-D)

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------

% set figure save subfolder
ddFigHelper.ResetParams()
figParams = ddFigHelper.GetAppData();
figParams.saveDir = fullfile(saveFiguresFolder,'LatchDecoder');
figParams.fig.fontSize = 8;
ddFigHelper.SetParams( figParams )

%% Figure 4(B): Example Trial
% plot gesture and latch decoder outputs to illustrate operation of the Latch decoder

sesI = 2; % choose session index
[sProbs,predState,latchState]  = CalculateProbabilities(sesData_GH_T11(sesI));

taskInfo = sesData_GH_T11(sesI).taskInfo;
selI = []; % indices of trials to plot (if empty, plot all)
trialNum = 9; % choose which trial to show in figure window
skipNoAction = true; % whether or not to show (false) or not show (true) no action state probabilities

PlotLatchProbs(sProbs,taskInfo, predState, latchState, selI, [], [],trialNum,skipNoAction)
ddFigHelper.SaveFigure(sprintf('Figure 4B - Toy example (T11), Session %d, Trial %d',sesI,trialNum))

% (note that PlotLatchProbs will plot all selI trials. With Matlab figure selected, use arrow keys (left/right) to browse trials)

% optionally add gesture legend to figure
% PlotGestureLegend(sProbs)
% ddFigHelper.SaveFigure(sprintf('Figure 4B - Toy example (T11), Session %d, Trial %d WITH LEGEND',sesI,trialNum))

%% Figure 4(C): Offline analysis of performance with Gesture only vs Latch decoder on Gesture Hero trials - T11
allXValPerSesT11(1,:) = ComparePerformanceAgainstLatchDecoder(sesData_GH_T11(1).feat,sesData_GH_T11(1).taskInfo); % session 1
allXValPerSesT11(2,:) = ComparePerformanceAgainstLatchDecoder(sesData_GH_T11(2).feat,sesData_GH_T11(2).taskInfo); % session 2

analysisTypeStr = 'postFirst'; % 'full', 'droppedPerSec, 'postFirst'
PlotXValGestureVsLatch(allXValPerSesT11,sesData_GH_T11(1),analysisTypeStr)

ddFigHelper.SaveFigure(sprintf('Figure 4C - Gesture Hero T11 - gesture vs latch per hold duration, %s',analysisTypeStr))

%% Figure 4(D): Offline analysis of performance with Gesture only vs Latch decoder on Gesture Hero trials - T5
allXValPerSesT5 = ComparePerformanceAgainstLatchDecoder(sesData_GH_T5.feat,sesData_GH_T5.taskInfo);

analysisTypeStr = 'postFirst'; % 'full', 'droppedPerSec, 'postFirst'
PlotXValGestureVsLatch(allXValPerSesT5,sesData_GH_T5,analysisTypeStr)

ddFigHelper.SaveFigure(sprintf('Figure 4D - Gesture Hero T5 - gesture vs latch per hold duration, %s',analysisTypeStr))





%% %%%%% HELPER FUNCTIONS %%%%% %%

function [sProbs,predState,latchState] = CalculateProbabilities(sesData)

    % Plot toy example
    feat = sesData.feat;
    taskInfo = sesData.taskInfo;

    %% Build latch decoder
    featInds = 1:384; % just TX and SP for this example trial
    trainOpt.performXVal = false;
    trainOpt.excludeTrialLen = [];
    trainOpt.taskParams.latch.noActionOffsetStartStop = [50 0];
    [mGesture, mLatch] = LatchedMultistateUpdated(taskInfo, feat(:,featInds),trainOpt);

    % hyperparams
    predOpt.gesturesThreshold = 0.9; %0.9998; % if empty use mean( prctile(gLiks,95) ); %0.85;
    predOpt.latchThreshold = 0.900;

    [predState,latchState, gLiks, lLiks] = PredictFromLatchDecoder(mGesture,mLatch,feat(:,featInds), predOpt);
    lLiks = cat(2,1-lLiks,lLiks); % Remove the no_action prob estimate

    %% Place into sProbs for plotting

    clear sProbs
    n = 1;
    sProbs(n).name = 'Attempt';
    sProbs(n).probs = lLiks;
    sProbs(n).labels = mLatch.info.train.uStates;
    sProbs(n).threshold = predOpt.latchThreshold;
    sProbs(n).colors = cat(1,zeros(1,3),ddHelper.decoder.colors(3,:)); %cat(1,zeros(1,3),cmap(length(mLatch.info.train.uStates)-1));

    n = n+1;
    sProbs(n).name = 'Gesture';
    sProbs(n).probs = gLiks;
    sProbs(n).labels = mGesture.info.train.uStates;
    sProbs(n).threshold = predOpt.gesturesThreshold;
    sProbs(n).colors = ddHelper.gestureColors; % cat(1,zeros(1,3),lines(length(mGesture.info.train.uStates)-1));
    
    if strcmp(sesData.participant,'t5') && strcmp(sesData.date,'2023.04.04')
        sProbs(n).colors = ddHelper.gestureColors([1 4 3 2 5 6 7 8],:);  % index finger down and power grasp are switchd in this session relative to T11 sessions
    end
    %%

    % tmp = CalcErrorSteps(taskInfo, gLiks, lLiks,feat(:,featInds), predOpt);

end

function PlotGestureLegend(sProbs)
    %% Legend (optional)
    gestIndex = find(ismember({sProbs.name},'Gesture'));
    
    fh = gcf;
    legLabels = sProbs(gestIndex).labels;
    if skipNoAction; legLabels(1) = []; end
    legInds = length(fh.Children(3).Children) - ((1:length(legLabels))*2 - 2);
    leg = legend([fh.Children(3).Children(legInds)],legLabels,'FontSize',4.5,'Location','best');

end




function allXVal = ComparePerformanceAgainstLatchDecoder(feat,taskInfo,featInds,gestures)

if nargin < 3
    featInds = 1:size(feat,2);%1:384;
end
if nargin < 4
    uGestLabels = unique(taskInfo.labels);
    gestures = uGestLabels;
end

uTrialDuration = [50 100, 200];

k = 5;

nIter = 10;

clear allXVal
predOpt.gesturesThreshold = 0.9; % if empty use mean( prctile(gLiks,95) ); %0.85;
predOpt.latchThreshold = 0.9;

allPerfLatch = {};
allPerfGest = {};
for di = 1:length(uTrialDuration)

    dur = uTrialDuration(di);
    dur_sec = dur./50;
    fprintf('\n%.0f seconds\n', dur_sec)

    isDur = ismember( taskInfo.trialDurationRounded, dur ); 
    isGest = ismember( taskInfo.labels,  gestures);  %select specific gestures if requested
    selI = find(isDur & isGest);
    selNoAction = selI;

    %% Train
    clear full postFirst
    lprintf()
    for xvi = 1:nIter
        lprintf('%02d of %02d\n', xvi,nIter)
        mv = mod(xvi-1,k);
        isReset = ~mv;
        xValNum = mv+1;
        if isReset
            xvalInds = crossvalind('KFold',length(selI),k);
        end
        trainInds = xvalInds ~= xValNum;
        testInds =  xvalInds == xValNum;


        trainOpt.excludeTrialLen = 300; % Steps
        trainOpt.selectAttempt = selI(trainInds); % Logical into startStops, labels
        trainOpt.selectNoAction = selNoAction(trainInds); % Logical into noActionStartStops
        trainOpt.performXVal = false;
        trainOpt.verbosity = 0;
        trainOpt.minMaxFeat = [0 400]; %make [] for no MRMR taking top feats
        trainOpt.sigAlpha = 0.001; % make =1 for no feature selection     
        [mGesture, mLatch] = LatchedMultistateUpdated(taskInfo, feat(:,featInds),trainOpt);

        testBufferWin = [-100 25];
        selTest = selI(testInds);
        testStartStops = taskInfo.startStops(selTest,:);
        testStartStopsWBuffer = testStartStops + testBufferWin;
        predWithBuffer = RowColon(testStartStopsWBuffer);
        testPredInds = RowColon(testStartStops);
        testPredFromBuffer = ismember(predWithBuffer,testPredInds); % May be issue if we wrap back into a prev trial

        [predState,latchState, gLiks, lLiks] = PredictFromLatchDecoder(mGesture,mLatch,feat(predWithBuffer,featInds), predOpt);
        predLabels = predState(testPredFromBuffer);


        %%
        [gv, gesturePred] = max(gLiks,[],2);
        gesturePred(gv<predOpt.gesturesThreshold) = 1;
        gesturePredLabels = gesturePred(testPredFromBuffer);
        %%

        testTask = GetTrainTest(taskInfo,selTest);

        % Update to relative startStops
        testLen = diff(testTask.startStops,[],2)+1;
        relStarts = cumsum([1; testLen]);
        relStartStops = cat(2, relStarts(1:end-1), relStarts(2:end)-1);
        testTask.startStops = relStartStops;
        [tmpL,prctLatch,nDroppedLatch] = CalcErrorSteps(testTask,predLabels);
        [tmpG,prctGest,nDroppedGesture] = CalcErrorSteps(testTask,gesturePredLabels);

        postFirst.latch(xvi) = 100-mean(prctLatch);
        postFirst.gesture(xvi) = 100-mean(prctGest);

        droppedPerSec.latch(xvi) = mean(nDroppedLatch)/dur_sec;
        droppedPerSec.gesture(xvi) = mean(nDroppedGesture)/dur_sec;
        %%

        testLabels = GetLabelsForStartStops(taskInfo.labels(selTest),testStartStops);
        testLabels = double(testLabels);

        isCorrect = (testLabels) == predLabels;
        xPrctLatch = sum(isCorrect)/length(isCorrect)*100;

        isCorrect = (testLabels) == gesturePredLabels;
        xPrctGesture = sum(isCorrect)/length(isCorrect)*100;

        full.latch(xvi) = xPrctLatch;
        full.gesture(xvi) = xPrctGesture;
        
        usedFeats.latch{xvi} = mLatch.params.featInds; %note that these are the used feats for attempt only decoder - not full latch
        usedFeats.gesture{xvi} = mGesture.params.featInds;

    end

    allXVal(di).postFirst = postFirst;
    allXVal(di).full = full;
    allXVal(di).droppedPerSec = droppedPerSec;
    allXVal(di).usedFeats = usedFeats;

end

end


function PlotXValGestureVsLatch(allXVal,sesData,typeStr,gestures)
% typeStr = 'postFirst'; % postFirst 'full'; droppedPerSec

if nargin < 4
    gestures = [];
end

figWH = [1.38 1.4];
ddFigHelper.CreateFigure([],figWH)
clf
axArgs = ddFigHelper.GetDefaultAxArgs();
axArgs.marg_h = [0.28 0.05];
axArgs.marg_w = [0.24 0.05];
ax = ttight_subplot(1,1,axArgs);

% uTrialDuration = unique(taskInfo.trialDurationRounded);
uTrialDuration = [50 100, 200];

plotVal = [];

nDecoders = 2;
clrs = ddHelper.decoder.colors;%lines(nDecoders);
ms = 20;

%%
ax = gca;
baseValue = 0;
oi = 0.2;
xs = linspace(-oi,oi,nDecoders);
for ii = 1:size(allXVal,2)
    v = [allXVal(:,ii).(typeStr)];
    [u(1),ci(1,:)] = GetConfidenceInterval([v.gesture],'ci'); % estimate CIs from bootstrapping 10-fold Xval results
    [u(2),ci(2,:)] = GetConfidenceInterval([v.latch],'ci');

    for jj = 1:nDecoders
        bx = ii + xs(jj);
        bh(jj) = bar(bx,u(jj),'BaseValue',baseValue);
        bh(jj).BarWidth = (1-oi-0.1)/2;
        bh(jj).FaceColor = clrs(jj,:);
    end

    plotCIs = false;
    if plotCIs
        for jj = 1:nDecoders
            errX = bh(jj).XEndPoints;% ii + xs(jj);    
            eb = errorbar(errX,u(jj),ci(jj,1),ci(jj,2),'.k', 'MarkerEdgeColor', 'k', 'Color','k');
            eb.LineWidth = 0.1;
        end
    end
end


for jj = 1:nDecoders
    hh(jj) = plot(nan,'sk','MarkerFaceColor',clrs(jj,:),'MarkerSize', ms);
end
xlim([0 length(allXVal)] + 0.5)
% bar(plotVal)
SetTickLabels(sprintfc('%.0f', uTrialDuration./50))
xlabel('Hold duration (sec)')
switch typeStr
    case {'full','postFirst'}
        ylabel('Percent Correct (%)')
    case 'droppedPerSec'
        ylabel('Drops/sec')
end
% lh = legend(hh,'Gesture','Latch');
% lh.Position(1:2) = [0.263 0.864];

% shrink yLim for better legibility in plotting
if ~strcmp(typeStr,'droppedPerSec')
    if strcmp(sesData(1).participant,'t11')
        if length(gestures)==3
            ylim([30 90])
        else
            ylim([30 80])
        end
    elseif strcmp(sesData(1).participant,'t5')
        ylim([50 100])
    end
end

end


function [trainTask, testTask] = GetTrainTest(taskInfo,selTrain,selTest)
if nargin < 3
    selTest = [];
end
nTrials = size(taskInfo.startStops,1);

if length(selTrain) ~= size(taskInfo.startStops,1)
    tmp = false(nTrials,1);
    tmp(selTrain) = true;
    selTrain = tmp;
end
if isempty(selTest)
    selTest = ~selTrain;
end

%%
trainTask = SelectTrialsFromTask(taskInfo,selTrain);
testTask  = SelectTrialsFromTask(taskInfo,selTest);

end


function taskInfo = SelectTrialsFromTask(taskInfo,selI)
    nTrials = length(selI); % size(taskInfo.startStops,1);
    fn = fieldnames(taskInfo);
    for ii = 1:length(fn)
        if size(taskInfo.(fn{ii}),1) == nTrials
            taskInfo.(fn{ii}) = taskInfo.(fn{ii})(selI,:);
        end
    end
end

function [numErrorAfterFirstDecode, prctError,nDroppedEvents] = CalcErrorSteps(taskInfo, gSweepProbs,lSweepProbs,feat, predOpt)
% [mv, predGestState] = max(gSweepProbs,[],2);
% predGestState(mv<0.95) = 1;
if nargin == 2
    predState = gSweepProbs;
else
    [predState,latchState] = PredictFromLatchDecoder(gSweepProbs,lSweepProbs,feat, predOpt);
end
%%
numErrorAfterFirstDecode = [];
prctError = [];
stateVars ={predState};% {predGestState,predState};
for trlI = 1:size(taskInfo.startStops,1)
    trlInds = RowColon(taskInfo.startStops(trlI,:));
    for jj = 1:length(stateVars)
        statVar = stateVars{jj};
    firstDecode = find(statVar(trlInds)>1,1);
    if isempty(firstDecode)
        firstDecode = 1; % Count all time steps as error if we never decoded
    end
    compareInds = firstDecode:length(trlInds);
    isCorrect = statVar(trlInds) == double(taskInfo.labels(trlI));
    numErrorAfterFirstDecode(trlI,jj) = sum(~isCorrect(compareInds)); % & latchState(trlInds(compareInds)));
    prctError(trlI,jj) = numErrorAfterFirstDecode(trlI,jj)/length(compareInds)*100;
    nDroppedEvents(trlI,jj) = sum( diff(~isCorrect)>0 );
    end
end
end


