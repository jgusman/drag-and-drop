%% Gesture Decoder (offline) vs. Latch Decoder (Fig 7)
% Performance metrics of Latch Decoder in Drag and Drop task with T11
% compared to simulated results of Gesture Decoder only

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------


% set figure save folder
ddFigHelper.ResetParams()
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'DD_GestureVsLatch'))

%% Get Drag and Drop Task Performance Metrics

dragTrialEpochs = {'Wait','Drag','Hold'};
out = GetDragAndDropPerformance(sesData_DD_T11,dragTrialEpochs);

%% FIGURE 7

plotFields = {'prctCorrectTrials', 'prctCorrectPostOnset', 'nDroppedEventsPerSec'};
plotType = 'Bar';
useMean = 1;


ddFigHelper.LogPrints('Drag and drop performance (mean).txt')

axArgs = ddFigHelper.GetDefaultAxArgs();
figWH = [6.48 1.8];
axArgs.gap = 0.13;
axArgs.marg_w = [0.11 0.01];
axArgs.marg_h = [0.12 0.15];

ddFigHelper.CreateFigure(382371,figWH)

clf
axs = ttight_subplot(1,length(plotFields),axArgs);

letterXYOffsets = [-0.11 -0.04];
for ii = 1:length(axs)
    plotField = plotFields{ii};
    axes(axs(ii))
    PlotPerformancePerEpoch(out, plotField,plotType, useMean, axs(ii))
    ddFigHelper.AddLetter(axs(ii),2+ii,letterXYOffsets(1),letterXYOffsets(2));
end

% legend (optional)
legH = gobjects();
legH(1) = bar(nan,'FaceColor',ddHelper.decoder.colors(1,:));
legH(2) = bar(nan,'FaceColor',ddHelper.decoder.colors(2,:));
lh = legend(legH,'Gesture decoder (offline)', 'Latch decoder','FontSize',7);
lh.Position(1:2) = [0.75    1-lh.Position(4)];

ddFigHelper.LogPrints()

ddFigHelper.SaveFigure('Figure 7 - Drag and Drop Performance (Mean, 95% bootci)')




%% Error Statistics

ddFigHelper.LogPrints('ErrorStatisticsGestureVsLatch.txt')

metric = 'nDroppedEvents'; %
[gestVals,~] = CombineAcrossSessions(out,metric);
gestOnly = [gestVals{1}{:}];
gestLatch = [gestVals{2}{:}];
totalErrorEvents(1) = sum(gestOnly,'all','omitnan');
totalErrorEvents(2) = sum(gestLatch,'all','omitnan');
fprintf('%s, Gest Only: %d \n', metric,totalErrorEvents(1))
fprintf('%s, Latch: %d \n', metric,totalErrorEvents(2))


metric = 'nErrorEvents'; %
[gestVals,~] = CombineAcrossSessions(out,metric);
gestOnly = [gestVals{1}{:}];
gestLatch = [gestVals{2}{:}];
totalErrorEvents(1) = sum(gestOnly,'all');
totalErrorEvents(2) = sum(gestLatch,'all');
fprintf('%s, Gest Only: %d \n', metric,totalErrorEvents(1))
fprintf('%s, Latch: %d \n', metric,totalErrorEvents(2))


metric = 'nErrorEventsLatched'; %
[gestVals,~] = CombineAcrossSessions(out,metric);
gestOnly = [gestVals{1}{:}];
gestLatch = [gestVals{2}{:}];
totalErrorEvents(1) = sum(gestOnly,'all');
totalErrorEvents(2) = sum(gestLatch,'all');
fprintf('%s, Gest Only: %d \n', metric,totalErrorEvents(1))
fprintf('%s, Latch: %d \n', metric,totalErrorEvents(2))

medianMetrics = {'errorDurationSteps','onset_steps2','offset_steps2'};

fprintf('\n\n')
for i = 1:length(medianMetrics)

    metric = medianMetrics{i}; %
    
    [gestVals,~] = CombineAcrossSessions(out,metric);
    gestOnly = [gestVals{1}{:}];
    gestLatch = [gestVals{2}{:}];
    if iscell(gestOnly(1))
        gestOnly = [gestOnly{:}];
        gestLatch = [gestLatch{:}];
    end
    
    % Wilcoxon Rank Sum test between simulated gesture only and actual latch decoder error metrics
    p = ranksum(gestOnly,gestLatch);

    avgErrorEvents(1) = median(gestOnly,'omitnan');
    avgErrorEvents(2) = median(gestLatch,'omitnan');
    
    fprintf('Median %s, Gest Only: %d \n', metric,avgErrorEvents(1))
    fprintf('Median %s, Latch: %d \n', metric,avgErrorEvents(2))
    fprintf('RS Test: p = %0.5f. \n', p)
end

ddFigHelper.LogPrints()




%% %%%%% HELPER FUNCTIONS %%%%% %%


function PlotPerformancePerEpoch(out,plotField,plotType,useMean,ax)
if ~exist('ax','var') || isempty(ax)
    figure(328221)
    clf
    axs = ttight_subplot(1,1);
    ax = axs(1);
end

[plotFieldStr, plotFields] = GetPlotFieldInfo();
plotStr = plotFieldStr.(plotField);

padDecoderStr = pad({'Gesture','Latch'});

epochPeriods = out(1).info.dragEpochs;


[gestVals,gestDelta] = CombineAcrossSessions(out,plotField);
%%

nInvalids = [];
clrs = ddHelper.decoder.colors;%lines(2);
ov = 0.2;
offsets = [-ov ov];


fprintf('\n%s\n', plotStr)
for axI = 1:length(epochPeriods)
    for decI = 1:2
    
        pd = gestVals{decI}{axI};
        isInvalid = isnan(pd);
        pd = pd(~isInvalid);
        nInvalids(axI,decI) = sum(isInvalid);
    
        if useMean
            u = mean(pd);
            try
                ci = bootci(1e3, {@mean, pd});
            catch
                ci = [u u];
            end
            ci = ci-u;
        else %(median)
            tmp = prctile(pd,[50 25 75]);
            u = tmp(1);
            ci = tmp(2:3);
            ci = ci-u;
        end    
    
        if useMean
            avgStr = 'mean';
            errorStr = '95 CI';
        else
            avgStr = 'median';
            errorStr = 'IQR';
        end
        fprintf('%s %s  %4.1f %s (%+02.1f %+02.1f %s)\n',epochPeriods{axI},padDecoderStr{decI}, u, avgStr, ci(1),ci(2), errorStr)
        
        errorX = axI + offsets(decI);
        c = clrs(decI,:);
        
        switch plotType
            case 'Bar'
                bh = bar(ax,errorX,u,'FaceColor',c);
                bh.BarWidth = ov*1.5;
                if length(pd) > 10
                    errorbar(ax,errorX,u,ci(1),ci(2),'.k','MarkerFaceColor',c,'MarkerEdgeColor','k');
                else
                    ptO = ov/2;
                    ptX = linspace(-ptO,ptO,length(pd)) + errorX;
                    plot(ax,ptX,pd,'.k','MarkerSize',5);
                    txtY = u; % max(pd);
                    text(errorX,txtY,sprintf('%.1f', u), 'HorizontalAlignment','center','VerticalAlignment','bottom');
                end
            case 'Error'
                % Only error bar
                errorbar(ax,errorX,u,ci(1),ci(2),'ok','MarkerFaceColor',c,'MarkerEdgeColor','k');
        end
    end

end
ax.YLabel.String = plotStr;
ax.XLim = [0.5 length(epochPeriods) + 0.5];
SetTickLabels(ax,epochPeriods)

end



function [plotFieldStr, plotFields] = GetPlotFieldInfo()

plotFieldStr.onset_steps = 'Decode reaction time';
plotFieldStr.nErrorEvents = '# incorrect gesture decode events';
plotFieldStr.nErrorSteps = '# incorrect gesture decode steps';
plotFieldStr.errorDurationSteps = 'error durations (steps)';
plotFieldStr.nErrorEventsLatched = '# of incorrect gesture decodes longer than 400ms';
plotFieldStr.nDroppedEvents = '# dropped events';
plotFieldStr.nDroppedEventsPerSec = '# drops/sec';
plotFieldStr.dropRateEvents = 'Expected drop rate (sec)';


plotFieldStr.droppedEventInds = 'Dropped events (sec)';
plotFieldStr.nCorrectSteps = '# correct steps';
plotFieldStr.prctCorrectSteps = 'Percent correct';
plotFieldStr.prctCorrectTrials = 'No error trials (%)';
plotFieldStr.prctCorrectPostOnset = 'Percent correct';
plotFieldStr.prctErrorSteps = 'Percent incorrect';
plotFieldStr.offset_steps = 'Offset time (steps)';
plotFieldStr.onset_steps2 = 'Onset time (steps)';
plotFieldStr.offset_steps2 = 'Offset time (steps)';


plotFields = fieldnames(plotFieldStr);
end


function [gestVals,gestDelta] = CombineAcrossSessions(out,plotFieldType, gestTypes)
if nargin < 3
    gestTypes = {'gestOnly','gestLatch'};
end
    nSes = length(out);
    
    %%
    gestVals = {};
    [~,plotFieldName] = GetFactorPerTrial(out,1,1,gestTypes{1},plotFieldType);
    for gtI = 1:length(gestTypes)
        gestType = gestTypes{gtI};
        nEpochs = size(out(1).(gestType).(plotFieldName),2);
        vals = cell(nEpochs,1);
        for ii = 1:nSes
            for jj = 1:nEpochs
                val = GetFactorPerTrial(out,ii,jj,gestType,plotFieldType);
                vals{jj} = cat(1, vals{jj}, val);
            end

        end

        gestVals{gtI} = vals;
    end
    
    nHist = 15;
    gestDelta = {};
    try
        for jj = 1:nEpochs
            [~,edges] = histcounts(cat(1,gestVals{1}{jj},gestVals{2}{jj}),nHist);
            n1 = histcounts(gestVals{1}{jj},edges);
            n2 = histcounts(gestVals{2}{jj},edges);
            n = n2-n1;
            edgeCenters = edges(1:end-1) + diff(edges)/2;
            tmp = arrayfun(@(x,y) repmat(x,y,1),edgeCenters,n,'UniformOutput',0);
            gestDelta{jj} = cat(1,tmp{:});
        end
    end
end


function [val,plotField] = GetFactorPerTrial(out,ii,jj,gestType, plotField)
trialLen_steps = out(ii).info.trialLength(:,jj);
trialLen_sec = trialLen_steps./50;
nTrials = length( trialLen_steps );

switch plotField
    case 'prctCorrectTrials'
        plotField = 'nDroppedEvents';
        % Convert to per block percent correct
        blks = out(ii).info.relCLBlock(:,jj);
        uBlks = unique( blks );
        for kk = 1:length(uBlks)
            selK = ismember( blks, uBlks(kk) );
            isDropped = out(ii).(gestType).(plotField)(selK,jj);
            isDropped(isnan(isDropped)) = 1; % Nan means always dropped. 
            isCorrect = ~logical(isDropped);
            val(kk,1) = sum(isCorrect)/length(isCorrect)*100;
        end
    case 'dropRateEvents'
        % Inverse of nDroppedEventsPerSec
        plotField = 'nDroppedEvents';
        factor = out(ii).(gestType).(plotField)(:,jj);
        factor(isnan(factor)) = 1; % Nan means always dropped. 
        val = trialLen_sec;
        update = factor~=0;
        val(update) = val(update)./factor(update);
    case 'nDroppedEventsPerSec'
        factor = 1./(trialLen_sec);
        plotField = 'nDroppedEvents';
        val = out(ii).(gestType).(plotField)(:,jj);
        val(isnan(val)) = 1; % Nan means always dropped. 
        val = val.*factor;
    case 'droppedEventInds'
        factor = repmat(1./50,nTrials,1);
        val = out(ii).(gestType).(plotField)(:,jj);
        val = arrayfun(@(x,y) cat(2,x{:}).*y, val,factor,'UniformOutput',0);
        val = cat(2, val{:})';
    case 'prctCorrectSteps'
        factor = 1./trialLen_steps*100;
        plotField = 'nCorrectSteps';
        val = out(ii).(gestType).(plotField)(:,jj);
        val = val.*factor;
    case 'prctErrorSteps'
        factor = 1./trialLen_steps*100;
        plotField = 'nErrorSteps';
        val = out(ii).(gestType).(plotField)(:,jj);
        val = val.*factor;
    otherwise
        val = out(ii).(gestType).(plotField)(:,jj);
end

end






function out = GetDragAndDropPerformance(sesData,dragEpochs)

% Note: In the DD task info structs here, each trial stage epoch (e.g., Wait, Drag, Hold, etc.) is a counted as a "trial" in the task

nSessions = length(sesData);
dragEpochsAll = {'Wait','Drag','Hold','Intertrial'};

if nargin < 2
    dragEpochs = dragEpochsAll;
end

for sesI = 1:nSessions

    taskInfo = sesData(sesI).taskInfo;
    outCL = sesData(sesI).outCL;

    decodedState = outCL.decodedState;   % decoded state at every 20ms time step in session
    gestureProbs = outCL.stateLL;        % gesture probabilities at every 20ms time step in session (from gesture decoder)
    gestureThresh = outCL.stateThresh;   % probability threshold at every 20ms time step in session (most likely same throughout)

    isOnTarget = outCL.isOnTarget;       % is cursor contacting cued target (checked every 20ms)

    uGest = unique( taskInfo.cuedGesture );
    no_action = uGest(1);
    selTrls = taskInfo.isWait & taskInfo.isDragAttempt & taskInfo.isGestureDecoderCL;
    selTrls = selTrls & ~ddHelper.ExcludePerformanceTrialsWithNS5Outliers(taskInfo); % remove trials excluded due to noise (NS5 outliers)
    dragEpochStartInds = find(selTrls);
    dragEpochRelInd = find(ismember(dragEpochs,dragEpochsAll)) - 1; %0:2

    uBlocks = unique(taskInfo.blockNumber(dragEpochStartInds));

    for esI = 1:length(dragEpochStartInds)
        
        label = taskInfo.cuedGesture(dragEpochStartInds(esI)); % gesture label for this trial

    
        %% Get decoder onsets and offsets (aka reaction times
        trlI = dragEpochStartInds(esI) + [-1 0 1 2 3 4]; %move, wait (and hold), drag, hold, intertrial, next trial's move period
        try
            trlInds = RowColon(taskInfo.startStops(trlI,:));
        catch
            trlI(end) = []; %last trial wont have a move state for the next trial
            trlInds = RowColon(taskInfo.startStops(trlI,:));
        end

        latchDecode = decodedState(trlInds); % decoded gesture state at each 20ms timesteps (what happened in CL using Latch decoder)

        % Get simulated Gesture decoder decodes (i.e. thresholded probabilities from Gesture decoder)
        blk = taskInfo.blockNumber(trlI(1));
        blkI = sesData(sesI).blocks==blk;

        [mv,mi] = max(gestureProbs(trlInds,:),[],2);                            % Find gesture state with highest likelihood at each timepoint of trial
        mi(mv<gestureThresh(trlInds)) = 1;                                      % Only count decode if above threshold (if not, default to 1 (no action))
        imagery = unique({taskInfo.actionMapPerBlock{blkI}.imagery},'stable');  % Get map numerical state to gesture name (categorical)
        gestOnlyDecode = categorical(mi,1:length(imagery),imagery);             % Change to categorical

        subfields = {'gestLatch', 'gestOnly'};
        decodeVars = {latchDecode, gestOnlyDecode};

        for decI = 1:length(subfields)  % for each decoder type (latch decoder and gesture only decoder)
            subfield = subfields{decI};

            ds = decodeVars{decI};

            [onset_steps2, offset_steps2] = GetGestureDecodeOnsetAndOffsetTimes(taskInfo,dragEpochStartInds, esI, trlI, trlInds, ds,isOnTarget, label);

            out(sesI).(subfield).onset_steps2(esI,1) = onset_steps2;
            out(sesI).(subfield).offset_steps2(esI,1) = offset_steps2;
        end


        %% Get other metrics by trial stage

        dragInds = taskInfo.startStops(dragEpochStartInds(esI),1):taskInfo.startStops(dragEpochStartInds(esI)+dragEpochRelInd(end),2);
        latchDecode = decodedState(dragInds);
        isDragging = latchDecode == label; % Why did this grow during the hold period?
        consecutiveStepsDragging = FindConsecutiveOnes(isDragging);

        for relI = 1:length(dragEpochRelInd)

            trlI = dragEpochStartInds(esI) + dragEpochRelInd(relI);
            trlInds = RowColon(taskInfo.startStops(trlI,:));

            % Should we look behind (early decodes?)
            latchDecode = decodedState(trlInds);
            epochLabel = taskInfo.cuedGesture(trlI);

            % Get Gesture state probability
            blk = taskInfo.blockNumber(trlI);
            blkI = sesData(sesI).blocks==blk;

            [mv,mi] = max(gestureProbs(trlInds,:),[],2);                            % Find highest gesture likelihood
            mi(mv<gestureThresh(trlInds)) = 1;                                      % Only count decode if above threshold (if not, default to 1 (no action)
            imagery = unique({taskInfo.actionMapPerBlock{blkI}.imagery},'stable');
            gestOnlyDecode = categorical(mi,1:length(imagery),imagery);             % Change to categorical

            subfields = {'gestLatch', 'gestOnly'};
            decodeVars = {latchDecode, gestOnlyDecode};
            for decI = 1:length(subfields)
                subfield = subfields{decI};

                [onset_steps, nErrorEvents, nErrorSteps, nDroppedEvents, droppedEventInds, nCorrectSteps, prctCorrectPostOnset, errorDurationSteps,nErrorEventsLatched,offset_steps] = GetEpochPerformanceMetrics(decodeVars{decI},epochLabel,relI, no_action);

                out(sesI).(subfield).onset_steps(esI,relI) = onset_steps;
                out(sesI).(subfield).nErrorEvents(esI,relI) = nErrorEvents;
                out(sesI).(subfield).nErrorSteps(esI,relI) = nErrorSteps;
                out(sesI).(subfield).errorDurationSteps{esI,relI} = errorDurationSteps;
                out(sesI).(subfield).nErrorEventsLatched(esI,relI) = nErrorEventsLatched;
                out(sesI).(subfield).nDroppedEvents(esI,relI) = nDroppedEvents;
                out(sesI).(subfield).droppedEventInds{esI,relI} = droppedEventInds;
                out(sesI).(subfield).nCorrectSteps(esI,relI) = nCorrectSteps;
                out(sesI).(subfield).prctCorrectPostOnset(esI,relI) = prctCorrectPostOnset;
                out(sesI).(subfield).offset_steps(esI,relI) = offset_steps;

            end

            out(sesI).info.sessionIndex(esI,relI) = sesI;
            out(sesI).info.taskEpochIndex(esI,relI) = trlI;
            out(sesI).info.trialLength(esI,relI) = length(trlInds);
            out(sesI).info.label(esI,relI) = epochLabel;
            out(sesI).info.relCLBlock(esI,relI) = find(uBlocks==blk);
            out(sesI).info.consecutiveStepsDragging{esI,relI} = consecutiveStepsDragging;

            out(sesI).info.dragEpochs = dragEpochs;

        end
    end
end

end


function [onset_steps, nErrorEvents, nErrorSteps, nDroppedEvents, droppedEventInds, nCorrectSteps, prctCorrectPostOnset,errorDurationSteps,nErrorEventsLatched,offset_steps] = GetEpochPerformanceMetrics(ds,label,relI, no_action)
% Onset
% ds = decoded state
% label = cued gesture
% relI = trial period index (e.g. Wait Drag Hold)
% no_action = the no action state / label

        onsetInd = find(ds==label,1);
        
        isWaitStage = relI == 1; % Only look at onset for wait stage
        if isempty(onsetInd) || ~isWaitStage
            onset_steps = nan;
        else
            onset_steps = onsetInd;
        end

        offsetInd = find(ds~=label,1);

        isIntertrialStage = relI == 4; % Only look at onset for wait stage
        if isempty(offsetInd) || ~isIntertrialStage || offsetInd == 1  %if first ind is error decode, then wont count it. (only count when intertrial period is entered with correct decode
            offset_steps = nan;
        else
            offset_steps = offsetInd;
        end
        
        % Errors
        % Number of times the wrong gesture is initially decoded
        isErrorDecode = ~ismember(ds,[no_action label]);
        nErrorSteps = sum(isErrorDecode);
        errorEventInds = find(diff([false; isErrorDecode])>0);
        errorEventEndInds = find(diff([isErrorDecode; false])<0);
        nErrorEvents = length( errorEventInds );
        
        errorDurationSteps = [];
        for i = 1:nErrorEvents
            errorDurationSteps(i) = length(errorEventInds(i):errorEventEndInds(i));
        end

        nErrorEventsLatched = sum(errorDurationSteps > 20); % i.e. greater than 400 ms

        nCorrectSteps = sum(ds == label);
        
        
        % Number of times dropped
        % For gesture + latch decoder
        % For gesture only
        if isempty(onsetInd)
            % Can't drop if you don't pick up.
            nDroppedEvents = nan;
            droppedEventInds = [];
            prctCorrectPostOnset = 0;
        else
            relInds = onsetInd:length(ds);
             
            isNotDragging = ~ismember(ds(relInds),label);
            if isWaitStage
                prctCorrectPostOnset = sum(~isNotDragging)/length(isNotDragging)*100;
            else
                prctCorrectPostOnset = sum(ismember(ds,label))/length(ds)*100;
            end
            
            nDroppedEvents = sum( diff(isNotDragging)>0 );
            
            droppedEventOnset = find([false; diff(isNotDragging)>0]);
            droppedEventInds = relInds(droppedEventOnset);
        end
        
end

function  [onset_steps2, offset_steps2] = GetGestureDecodeOnsetAndOffsetTimes(taskInfo, dragEpochStartInds, esI, trlI,trlInds, ds, isOnTarget, label)

    % find data indices of of Wait, Intertrial, and Drag epochs relative to the start of the Move epoch
    relWaitStartInd = taskInfo.startStops(dragEpochStartInds(esI),1) - taskInfo.startStops(dragEpochStartInds(esI)-1,1); % beginning of wait period (relative to start of move to target trial stage)
    relIntertrialStartInd = taskInfo.startStops(dragEpochStartInds(esI)+3,1) - taskInfo.startStops(dragEpochStartInds(esI)-1,1);  % beginning of intertrial period (relative to start of move to target trial stage)
    relIntertrialEndInd = taskInfo.startStops(dragEpochStartInds(esI)+3,2) - taskInfo.startStops(dragEpochStartInds(esI)-1,1);  % end of intertrial period (relative to start of move to target trial stage)
    relDragStartInd = taskInfo.startStops(dragEpochStartInds(esI)+1,1) - taskInfo.startStops(dragEpochStartInds(esI)-1,1); % beginning of wait period (relative to start of move to target trial stage)
    
    relIsOnTarget = isOnTarget(trlInds); % onTarget index relative to start of move to target trial stage
    onTargetEvents = find(diff(relIsOnTarget)>0)+1;
    
    onTargetEventInd = find(onTargetEvents<=relWaitStartInd,1,'last'); % get last on target event before beginning of Wait period
    onTargetInd = onTargetEvents(onTargetEventInd);

    % look for gesture decoding onset relative to the last On Target index (not the beginning of the wait period!)
    if ds(onTargetInd) ~= label  %if not already doing the gesture once reaching the target (wait stage)
        indexInd = onTargetInd;
        while indexInd <= relDragStartInd && ds(indexInd) ~= label
            indexInd = indexInd + 1;
        end
        onset_steps2 = indexInd - onTargetInd;
        if indexInd >= relDragStartInd
            onset_steps2 = NaN; % if onset doesnt occur during wait period
        end
    elseif ds(onTargetInd) == label  %if already doing the gesture when reaching the target (wait stage)
        indexInd = onTargetInd;
        while indexInd > 0 && ds(indexInd) == label
            indexInd = indexInd - 1;
        end
        onset_steps2 = indexInd - onTargetInd;  % will be negative number
    end
    
    % look for gesture decoding offset relative to the start of the intertrial period
    if ds(relIntertrialStartInd) == label  %if still doing the gesture once reaching the end of the hold stage (start of intertrial)
        indexInd = relIntertrialStartInd;
        while indexInd <= relIntertrialEndInd+50 && indexInd <= length(ds) && ds(indexInd) == label
            indexInd = indexInd + 1;
        end
        offset_steps2 = indexInd - relIntertrialStartInd;
        if indexInd >= relIntertrialEndInd+50
            offset_steps2 = NaN; % if offset doesnt occur during intertrial period
        end
    elseif ds(relIntertrialStartInd) ~= label  %if not doing the gesture once hit the intertrial period, look back to see when last time they were doing it
        indexInd = relIntertrialStartInd;
        while indexInd > relWaitStartInd && ds(indexInd) ~= label
            indexInd = indexInd - 1;
        end
        offset_steps2 = indexInd - relIntertrialStartInd;  % will be negative number
    end
    
    if diff(taskInfo.startStops(trlI(1),:)) >= 25*50-1  %if previous movement to target was a timeout (up to 25 sec long)
        onset_steps2 = NaN;
    end
    
    if diff(taskInfo.startStops(trlI(3),:)) >= 25*50-1  %if previous movement to target (drag epoch) was a timeout (up to 25 sec long)
        offset_steps2 = NaN;
    end

end
