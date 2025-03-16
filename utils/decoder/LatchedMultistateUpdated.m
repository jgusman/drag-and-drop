function [mGesture, mLatch, holdout] = LatchedMultistateUpdated(taskInfo, features, trainOpt)
%% Get task info out
%
% Should no action be split into prep/move/rest?
% 
%% Default Values
if nargin < 3
    trainOpt = [];
end
dfTrainOpt.verbosity = 1;
dfTrainOpt.taskToTrainFunction = @GetTrainInfoFromTaskInfo;
dfTrainOpt.excludeFeatInds = [];
dfTrainOpt.holdoutTrial_prct = 0; % Where 100 is 100%
dfTrainOpt.excludeOutlierThresh_prct = 5; 
dfTrainOpt.excludeMultitargetOnset = 2;
dfTrainOpt.excludeTrialLen = 400; % Steps
dfTrainOpt.selectAttempt = []; % Logical into startStops, labels
dfTrainOpt.selectNoAction = []; % Logical into noActionStartStops
dfTrainOpt.holdoutBlocks = [];
dfTrainOpt.trainBlocks = [];
dfTrainOpt.plotHoldAssessment = false;
dfTrainOpt.performXVal = true;

dfTrainOpt.sigAlpha = 0.001; 
dfTrainOpt.minMaxFeat = [];

dfTrainOpt.decoderParams.ldaRegularization = 0.001;
dfTrainOpt.decoderParams.dimReductionMethod = 'lda';
dfTrainOpt.decoderParams.smoothLikAmount = [15 0];
dfTrainOpt.decoderParams.xValVerbosity = 0;

dfTrainOpt.decoderParams.gesture.smoothDataAmount = [5 0];
dfTrainOpt.decoderParams.gesture.enforceTransitionMatrix = [];

dfTrainOpt.decoderParams.latch.smoothDataAmount = [15 0];
dfTrainOpt.decoderParams.latch.kFoldIter = 5;
dfTrainOpt.decoderParams.latch.enforceTransitionMatrix = 1-10^-13;


dfTrainOpt.taskParams.gesture.labelType      = 'Gesture';
dfTrainOpt.taskParams.gesture.actionOffsetStartStop = [0 0];
dfTrainOpt.taskParams.gesture.noActionOffsetStartStop = [0 0];

dfTrainOpt.taskParams.latch.labelType      = 'Latch';
dfTrainOpt.taskParams.latch.actionOffsetStartStop = [0 0];
dfTrainOpt.taskParams.latch.noActionOffsetStartStop = [0 0];
dfTrainOpt.taskParams.latch.noActionVal = 'no_action';

trainOpt = MergeBintoA(dfTrainOpt,trainOpt);


%% Convert taskInfo to trainInfo and get holdout task info
% [taskInfo_train, holdout] = GetTrainInfoFromTaskInfo(taskInfo, trainOpt);
[taskInfo_train, holdout] = GetTrainInfoFromGestureHeroTaskInfo(taskInfo, trainOpt);
%% Gesture decoder
% Set Up Gesture Task Info 
% Merge gesture specific task params into task

mGesture = TrainGestureDecoder(features, taskInfo_train, trainOpt);

%% Latch Decoder
mLatch = TrainLatchDecoder(features, taskInfo_train, trainOpt);

%% Hold out data

% Now train both decoder on all data
% mLatch.Train();
% mGesture.Train();

% Test on hold out data
% TestGestureLatchDecoder(features, holdout, mGesture, mLatch); % @Tommy is this obsolete?


%% X val for testing...
% XValDecoder(features, decoderParams)
end

function mGesture = TrainGestureDecoder(features, taskInfo_train, trainOpt)
    
    % Get Gesture-Specific Decoder Params
    gDecoderParams = MergeBintoA(trainOpt.decoderParams,trainOpt.decoderParams.gesture);
    
    
    gTaskInfo = MergeBintoA(taskInfo_train, trainOpt.taskParams.gesture);
    %
    if any(isundefined( gTaskInfo.labels ) )
        error('Found undefined labels!')
    end
    LocalPrintF(trainOpt,'Gesture decoder train labels (per trial)\n');
    if trainOpt.verbosity
        summary( gTaskInfo.labels )
    end
    
    LocalPrintF(trainOpt,'Updating start stop epochs\n')
    gTaskInfo = UpdateStartStops( gTaskInfo );
    
    if ~isfield(gTaskInfo, 'rankFeat')
        LocalPrintF(trainOpt,'Ranking features...\n')
        gTaskInfo.rankFeat = RankSigFeatures(features, gTaskInfo.labels, gTaskInfo.startStops, trainOpt.excludeFeatInds, trainOpt.sigAlpha, trainOpt.minMaxFeat );
    end
    LocalPrintF(trainOpt,'Done!\n')

    %% Train the decoder
    mGesture = TrainStateDecoder(features, gTaskInfo, gDecoderParams, trainOpt.performXVal);

    if trainOpt.performXVal
        LocalPrintF(trainOpt,'Decoder latch performance stats on last x-val fold.\n')
        AssessHoldError(mGesture, trainOpt.plotHoldAssessment);
    end

end


function mLatch = TrainLatchDecoder(features, taskInfo_train, trainOpt)
    
    % Get Latch-Specific Decoder Params
    lDecoderParams = MergeBintoA(trainOpt.decoderParams,trainOpt.decoderParams.latch);

    %% Specify Extra Task Info
    lTaskInfo = MergeBintoA(taskInfo_train,trainOpt.taskParams.latch);
    % lTaskInfo = ConvertLatchToOnsetOffset(lTaskInfo)
    %%
    % Change all labels to latch
    lTaskInfo.labels(~ismember(lTaskInfo.labels, lTaskInfo.noActionVal)) = 'Latch';
    lTaskInfo.labels = removecats(lTaskInfo.labels);
    %%
    
    %
    if any(isundefined( lTaskInfo.labels ) )
        error('Found undefined labels!')
    end
    LocalPrintF(trainOpt,'Latch decoder train labels (per trial)\n');
    if trainOpt.verbosity
        summary( lTaskInfo.labels )
    end
    
    LocalPrintF(trainOpt,'Updating start stop epochs\n')
    lTaskInfo = UpdateStartStops( lTaskInfo );
    if ~isfield(lTaskInfo, 'rankFeat')
        LocalPrintF(trainOpt,'Ranking features...\n')
        lTaskInfo.rankFeat = RankSigFeatures(features, lTaskInfo.labels, lTaskInfo.startStops, ...
                                    trainOpt.excludeFeatInds,trainOpt.sigAlpha, trainOpt.minMaxFeat);
    end
    LocalPrintF(trainOpt,'Done!\n')
    
    %% This works! Latch MCC: 0.634 ±0.014    Prct correct: 95.225%
    mLatch = TrainStateDecoder(features, lTaskInfo, lDecoderParams, trainOpt.performXVal);

    if trainOpt.performXVal
        LocalPrintF(trainOpt,'Decoder latch performance stats on last x-val fold.\n')
        AssessHoldError(mLatch, trainOpt.plotHoldAssessment);
    end
end

function lTaskInfo = ConvertLatchToOnsetOffset(lTaskInfo)
    isAction = lTaskInfo.isAction;
    startStops = lTaskInfo.startStops(isAction,:);
    onsetWin = [-25 25];
    offsetWin = [-25 25];
    noActionWin = [offsetWin(2) onsetWin(1)] + [1 -1];
    onsetStartStops = startStops(:,1) + onsetWin;
    offsetStartStops = startStops(:,2) + offsetWin;
    noActionStartStops = lTaskInfo.startStops(~isAction,:) + noActionWin;
    
    %noActionStartStops = [offsetStartStops(1:end-1,2)+1 onsetStartStops(2:end,1)-1];
%     noActionLen = diff(noActionStartStops,[],2);
%     figure(33); cla; plot(noActionLen)
%     rmvNoAction = noActionLen > 50*5;
    isOverlap = any(ismember( RowColon(noActionStartStops), cat(1,onsetStartStops,offsetStartStops)),2);
    if any(isOverlap)
        error('Overlapped attempt onset/offset with no action')
    end
    
    stStps = {noActionStartStops, onsetStartStops, offsetStartStops};
    lbls = {'no_action', 'Onset', 'Offset'};
    startStops = [];
    labels = [];
    for ii = 1:length(stStps)
        startStops = cat(1,startStops,stStps{ii});
        labels = cat(1,labels, repmat(ii, size(stStps{ii},1),1));
    end
    labels = categorical(labels,1:length(stStps),lbls);
    [~,si] = sort(startStops(:,1));
    startStops = startStops(si,:);
    labels = labels(si);
    
    lTaskInfo.startStops = startStops;
    lTaskInfo.labels = labels;
end


function XValDecoder(features, decoderParams)
    % I think this is for assessing different feats?
    %% Xval
    foldVals = 10:20:400; % linspace(0,0.9,15);
    % foldVals = 0:5:25;
    mccPerXval = [];
    withinWinCorrect = [];
    for xi = 1:length(foldVals)
        %%
    %     if foldVals(xi) == 0
    %         decoderP.enforceTransitionMatrix = [];
    %         transDiag = 1;
    %     else
    %         transDiag = 1 - 10^-xi;
    %         decoderP.enforceTransitionMatrix = transDiag;
    %     end
    transDiag = 1 - 10^-6; % 10^-xi;
    decoderP.enforceTransitionMatrix = transDiag;
    % latchDecoderP.maxComponents = slc.sSLC.decoders.multistate.maxNumPCs;
    decoderP.ldaRegularization = 0.1; % foldVals(xi);
    decoderP.smoothDataAmount = [15 0]; % [foldVals(xi) 0]; % [5, 0];
    decoderP.xValVerbosity = 0;
    decoderP.kFoldIter = 5;
    decoderParams.params.decoder.numFeat = foldVals(xi);
    decoderParams.selectedFeatInds = selectedFeatInds(1:foldVals(xi));
    
    mDecoder = TrainStateDecoder(features, decoderParams, 1, opt );
    
    mccPerXval(xi,:) = [mDecoder.performance.mcc, mDecoder.performance.std.mcc];
    
    holdErrOut = AssessHoldError(mDecoder,0);
    fn = fieldnames(holdErrOut);
    for fnI = 1:length(fn)
        withinWinCorrect(xi,fnI) = mean(holdErrOut.(fn{fnI}).withinCorrect_prct); 
    end
    
    
    end
    %%
    fprintf('\n\nCross validation results:\n')
    fprintf('Within window correct class order:  ')
    fprintf('%s | ', fn{:})
    fprintf('\n\n')
    %%
    for xi = 1:length(foldVals)
        withinWinCorrectStr = sprintf('%.2f%%, ', withinWinCorrect(xi,:));
        fprintf('%s xVal value: %.3f \t  %.3f ±%.3f MCC | Within window correct: %.3f%%\n', ...
            decoderParams.labelType, foldVals(xi), mDecoder.performance.mcc, mDecoder.performance.std.mcc, withinWinCorrectStr);
    
    end

end

function [trainInfo, holdout] = GetTrainInfoFromGestureHeroTaskInfo(taskInfo, opt)
    
    
    holdout = [];
    
    %% Selection
    if ~isempty(opt.selectAttempt)
        taskInfo.startStops = taskInfo.startStops(opt.selectAttempt,:);
        taskInfo.labels = taskInfo.labels(opt.selectAttempt,:);
    end
    
    if ~isempty(opt.selectNoAction)
        taskInfo.noActionStartStops = taskInfo.noActionStartStops(opt.selectNoAction,:);
    end
    
    %% Combine no action trials
    nNoActionTrials = size(taskInfo.noActionStartStops,1);
    nActionTrials = size(taskInfo.startStops,1);
    
    startStops = cat(1, taskInfo.startStops, taskInfo.noActionStartStops);
    labels = taskInfo.labels;
    labelNoActionInds = nActionTrials + (1:nNoActionTrials);
    labels(labelNoActionInds) = 'no_action'; %repmat('no_action', nNoActionTrials, 1)
    isAction = cat(1,true(nActionTrials,1), false(nNoActionTrials,1));
    
    [~,si] = sort(startStops(:,1));
    
    trainInfo.startStops = startStops(si,:);
    trainInfo.labels = labels(si);
    trainInfo.isAction = isAction(si);
    
    
    %% Exclusion criteria
    % Compute average based on ns5 outliers per start stop
    LocalPrintF(opt,'\n\n--Outlier and hold out summary--\n')
    avgOutlier = [];
    for ii = 1:length(trainInfo.startStops)
        avgOutlier(ii) = mean( taskInfo.prctNS5Outliers( RowColon(trainInfo.startStops(ii,:)) ) );
    end
    excludeOutliers = avgOutlier > opt.excludeOutlierThresh_prct;
    LocalPrintF(opt,'Found %d outliers out of %d trials\n', sum(excludeOutliers), length(excludeOutliers));
    excludeTrials = excludeOutliers(:);
    
    if ~isempty(opt.excludeTrialLen)
        trialLen = diff(trainInfo.startStops,[],2)+1;
        
        excludeLongTrials = trialLen > opt.excludeTrialLen;
        LocalPrintF(opt,'Found %d long trials (> %d)\n', sum(excludeLongTrials),opt.excludeTrialLen)
        excludeTrials = excludeTrials | excludeLongTrials(:);
    end
    
    
    trainInfo.startStops(excludeTrials,:) = [];
    trainInfo.labels(excludeTrials) = [];
    trainInfo.isAction(excludeTrials) = [];
end

function [trainInfo, holdout] = GetTrainInfoFromTaskInfo(taskInfo, opt)
    
    % N = 1:20;
    % table(taskInfo.trialStage(N), taskInfo.cuedGesture(N), taskInfo.isDragAttempt(N), taskInfo.verboseTrialStage(N))
    
    
    %% Exclusion criteria
    % Compute average based on ns5 outliers per start stop
    LocalPrintF(opt,'\n\n--Outlier and hold out summary--\n')
    avgOutlier = [];
    for ii = 1:length(taskInfo.startStops)
        avgOutlier(ii) = mean( taskInfo.prctNS5Outliers( RowColon(taskInfo.startStops(ii,:)) ) );
    end
    excludeOutliers = avgOutlier > opt.excludeOutlierThresh_prct;
    LocalPrintF(opt,'Found %d outliers out of %d trials\n', sum(excludeOutliers), length(excludeOutliers));
    excludeTrials = excludeOutliers(:);
    
    if opt.excludeMultitargetOnset
        % Exclude trials with multiple target entires, 
        % There may be poor control. 
        [nOnTargPerTrial_CO,nOnTargPerTrial_WT] = GetTargetOnsetPerTrial(taskInfo);
        excludeMultiTargetOnset = nOnTargPerTrial_CO > opt.excludeMultitargetOnset;
        excludeMultiTargetOnset = excludeMultiTargetOnset | nOnTargPerTrial_WT > opt.excludeMultitargetOnset;
    
        LocalPrintF(opt,'Excluding %d trials due to more than %d multi-target onsets\n', sum(excludeMultiTargetOnset), opt.excludeMultitargetOnset);
        excludeTrials = excludeTrials | excludeMultiTargetOnset(:);
    end
    
    excludeUndefined = isundefined(taskInfo.cuedGesture);
    excludeTrials = excludeTrials | excludeUndefined(:);
    
    LocalPrintF(opt,'Excluding %d trials in total\n', sum(excludeTrials));
    %%
    
    % Do not include drag trials or wait trials
    noActionTrialLogical = ~taskInfo.isDragAttempt & ~taskInfo.isWait & ~excludeTrials;
    actionTrialLogical = ((taskInfo.isGestureOnset & ~taskInfo.isClickAttempt) | taskInfo.isDragAttempt) & ~excludeTrials;
    
    gestureStartStops = taskInfo.startStops;
    noActionStartStops = taskInfo.startStops; 
    
    % Exclude the on-target during center out. As that is when gesture
    % click/drags start.
    centerOutTrials = find(taskInfo.isCenterOut);
    for ii = 1:length(centerOutTrials)
        coTrl = centerOutTrials(ii);
        if ~isempty(taskInfo.onTargetStartStops{coTrl})
            % Set the no-action center out trial stop to be just before the
            % first target onset.
            noActionStartStops(coTrl,2) = taskInfo.onTargetStartStops{coTrl}(1,1)-1;
            
            % For Attempts, 
            % Set the center-out trial-start to be just after the last target onset.
            gestureStartStops(coTrl,1) = taskInfo.onTargetStartStops{coTrl}(end,2)+1;
        else
            % Never made it on target. Exclude for the gesture logical
            actionTrialLogical(coTrl) = false;
        end
    end
    
    % Pull out the noAction trials
    noActionStartStops = noActionStartStops(noActionTrialLogical,:);
    
    % Pull out the attempt start stops
    gestureStartStops = gestureStartStops(actionTrialLogical,:);
    
    % Sanity check for overlap
    if ~isempty(intersect(RowColon(noActionStartStops), RowColon(gestureStartStops)))
        vprintf(2,'Overlapping Gesture and No action inds!\n\nStopping!! (keyboard), press quit debugging to exit debugging!\n')
        keyboard;
    end
        
    %
    maxTrialLen = 400;
    for ii = 1:size(gestureStartStops,1)
        if diff(gestureStartStops(ii,:)) > maxTrialLen
            gestureStartStops(ii,2) = gestureStartStops(ii,1) + maxTrialLen-1;
        end
    end
    for ii = 1:size(noActionStartStops,1)
        if diff(noActionStartStops(ii,:)) > maxTrialLen
            noActionStartStops(ii,2) = noActionStartStops(ii,1) + maxTrialLen-1;
        end
    end
    
    
    %% Merge trials...
    % aa = find( (noActionStartStops(2:end,1) - noActionStartStops(1:end-1,2)) == 1 )
    
    %% Combine Action and No Action
    isAction = cat(1, true(size(gestureStartStops,1),1), false(size(noActionStartStops,1),1));
    startStops = cat(1,gestureStartStops, noActionStartStops);
    labels = cat(1, taskInfo.cuedGesture(actionTrialLogical), taskInfo.cuedGesture(noActionTrialLogical));
    [~,si] = sort( startStops(:,1) );
    labels = labels(si);
    startStops = startStops(si,:);
    isAction = isAction(si);
    
    %% Exclude and hold out trials 
    
    
    
    
    
    % Combine all data into trainInfo
    trainInfo.startStops = startStops;
    trainInfo.labels = labels;
    trainInfo.labels(~isAction) = 'no_action';
    trainInfo.isAction = isAction;
    
    % Do we select train events based on blocks?
    selectTrainWithBlockNum = ~isempty(opt.trainBlocks);
    selectHoldoutWithBlockNum = ~isempty(opt.holdoutBlocks);
    if selectTrainWithBlockNum && selectHoldoutWithBlockNum
        commonBlocks = intersect(opt.trainBlocks,opt.holdoutBlocks);
        if ~isempty(commonBlocks)
            error(sprintf('Train and holdout block selection cannot have common blocks!\nCommon blocks: %s', commonBlocks))
        end
    end
    
    
    % Trial selection
    nOrgTrials = length(trainInfo.isAction);
    selTrainTrialInds = 1:nOrgTrials; % select all train trials
    
    if selectTrainWithBlockNum
%         selTrainTrialInds = find(all(ismember(taskInfo.blockNumber( GetGestureHoldEpochInds(taskInfo, trainInfo.startStops) ), opt.trainBlocks),2));
        selTrainTrialInds = find(all(ismember(taskInfo.blockNumber( trainInfo.startStops ), opt.trainBlocks),2));
    end
    
    if selectHoldoutWithBlockNum
%         holdoutTrialInds = find(all(ismember(taskInfo.blockNumber( GetGestureHoldEpochInds(taskInfo,trainInfo.startStops) ), opt.holdoutBlocks),2));
        holdoutTrialInds = find(all(ismember(taskInfo.blockNumber( trainInfo.startStops ), opt.holdoutBlocks),2));
        nHoldoutTrials = length(holdoutTrialInds);
    else
        nHoldoutTrials = round( nOrgTrials * (opt.holdoutTrial_prct/100) );
        holdoutTrialInds = (nOrgTrials - nHoldoutTrials + 1):nOrgTrials;
    end
    
    selTrainTrialInds = setdiff(selTrainTrialInds,holdoutTrialInds);
    nTrainTrials = length(selTrainTrialInds);
    
    
    holdout.isAction = trainInfo.isAction(holdoutTrialInds,:);
    holdout.labels = trainInfo.labels(holdoutTrialInds);
    holdout.startStops = trainInfo.startStops(holdoutTrialInds,:);
    
    % Remove non-selected trials
    
    trainInfo.isAction = trainInfo.isAction(selTrainTrialInds,:);
    trainInfo.labels = trainInfo.labels(selTrainTrialInds);
    trainInfo.startStops = trainInfo.startStops(selTrainTrialInds,:);
    
    %%
    % Print sumary
    LocalPrintF(opt,'\n--Input task info summary--\n')
    LocalPrintF(opt,'Gesture summary (steps)\n')
    % summary( taskInfo.stateCue( RowColon(gestureStartStops) ) )
    if opt.verbosity
        summary( trainInfo.labels( trainInfo.isAction )); %ismember( taskInfo.epochStartStops(:,1), (gestureStartStops)) ) )
    end
    % PrintNewLines(2)
    LocalPrintF(opt,'\n\n')
    LocalPrintF(opt,'No action summary (steps)\n')
    % summary( taskInfo.stateCue( RowColon(noActionStartStops) ) )
    if opt.verbosity
        summary( trainInfo.labels( ~trainInfo.isAction ));
    end
    
    
    LocalPrintF(opt,'\n\nOriginal valid trials: %d\n', nOrgTrials)
    LocalPrintF(opt,'Selected train trials: %d (%.0f%%)\n', nTrainTrials, (nTrainTrials/nOrgTrials)*100)
    LocalPrintF(opt,'Hold out trials:       %d (%.0f%%)\n', nHoldoutTrials, (nHoldoutTrials/nOrgTrials)*100);
    LocalPrintF(opt,'Hold out gesture/no-action trial summary\n')
    if opt.verbosity
        summary(holdout.labels)
    end
    LocalPrintF(opt,'\n')

end

%%

function [nOnTargPerTrial_CO,nOnTargPerTrial_WT] = GetTargetOnsetPerTrial(taskInfo)
    nOnTargs_CO = cellfun(@(x) size(x,1), taskInfo.onTargetStartStops(taskInfo.isCenterOut));
    nOnTargs_WT = cellfun(@(x) size(x,1), taskInfo.onTargetStartStops(taskInfo.isWait));
    
    nOnTargPerTrial_CO = nan(size(taskInfo.onTargetStartStops));
    nOnTargPerTrial_WT = nan(size(taskInfo.onTargetStartStops));
    for ii = 1:length(nOnTargs_CO)
        selBlkI = find(ismember(taskInfo.trialSetInd,ii));
        nOnTargPerTrial_CO(selBlkI) = nOnTargs_CO(ii);
        nOnTargPerTrial_WT(selBlkI) = nOnTargs_WT(ii);
    end
end


function decoderTaskInfo = UpdateStartStops( decoderTaskInfo )
    % Update start stops by offsets
    %   .params.task.actionOffsetStartStop   = [startOffset, stopOffset]
    %   .params.task.noActionOffsetStartStop = [startOffset, stopOffset]
    % 
    % Assumptions: 
    %   We do not have a fixed window size
    %   Labels per trial (start-stop) will be fixed.
    
    %% Update offsets
    isAction = decoderTaskInfo.isAction;
    decoderTaskInfo.startStops(isAction,:)  = decoderTaskInfo.startStops(isAction,:)  + decoderTaskInfo.actionOffsetStartStop;
    decoderTaskInfo.startStops(~isAction,:) = decoderTaskInfo.startStops(~isAction,:) + decoderTaskInfo.noActionOffsetStartStop;
    
    
    % Error if any trials that are <= 1
    isTrialLessThanOne = diff(decoderTaskInfo.startStops,[],2) <= 1;
    % if any(isTrialLessThanOne)
    %     error('Found %d trials less than 1', sum(isTrialLessThanOne))
    % end
    decoderTaskInfo.startStops( isTrialLessThanOne, : ) = [];
    decoderTaskInfo.labels( isTrialLessThanOne ) = [];
    decoderTaskInfo.isAction( isTrialLessThanOne ) = [];

end



function DebugPlotLabels(taskInfo)
    labels = double(taskInfo.labels);
    startStops = taskInfo.startStops;
    figure(2);
    cla
    plot(labels)
    plot(startStops(:,1), labels(startStops(:,1)), 'go')
    plot(startStops(:,2), labels(startStops(:,2)), 'rx')
    
    pltVal = nan(size(labels));
    pltVal(RowColon(noActionStartStop)) = 1;
    plot(pltVal, 'LineWidth', 3 )
    
    pltVal = nan(size(labels));
    pltVal(RowColon(startStops)) = labels(RowColon(startStops));
    plot(pltVal, 'LineWidth', 3 )

end

function starStops = GetOnTargetStartStops(taskInfo, isLast)
if nargin < 3
    isLast = 1;
end
onTargStartStops = taskInfo.onTargetStartStops;
starStops = [];
nOnTargs = [];
for ii = 1:length(onTargStartStops)
    if ~isempty(onTargStartStops{ii})
        nOnTargs(end+1) = size(onTargStartStops{ii},1);
        if isLast
            starStops(ii,:) = onTargStartStops{ii}(end,:); % Last on target onset
        else
            starStops(ii,:) = onTargStartStops{ii}(1,:); % First on target onset
        end
    else
        % Target onset is typically at the end of a trial
        starStops(ii,:) = [taskInfo.startStops(ii,2) taskInfo.startStops(ii,2)];
    end
end
  
end

function LocalPrintF(opt,varargin)
if opt.verbosity
    fprintf(varargin{:});
end
end
