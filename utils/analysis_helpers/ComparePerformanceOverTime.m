
function swp = ComparePerformanceOverTime(sesData,sweepInfo, gestures)
% swp = ComparePerformanceOverTime(sesData,sweepInfo)
% Cross validate LDA decoder at different offsets of a trial.
% 
% Sweep across attempt durations
% Sweep across trial offsets.
% 
% sweepInfo structure format
% 
% sweepInfo.attemptWindow = [-24 0]; % What goes into the lda for decoding. [-24 0] is 0.5 seconds back.
% sweepInfo.attemptOffset = [-50 100]; % Sweep boundaries relative to attempt trial [start stop]
% sweepInfo.noAttemptWindow = [25 150]; % No-attempt period relative to noActionStartStops start. (attempt trial end)
% sweepInfo.stepSize = 15;
% sweepInfo.featInds = 1:384;
% sweepInfo.enforceTransitionMatrix = 0; % 0 = No HMM. [] = Yes HMM. Transition matrix from train data trial transitions.


if nargin < 3
    uGestLabels = unique(sesData.taskInfo.labels);
    gestures = uGestLabels;
end


%% Get data and task info
    feat     = sesData.feat;
    taskInfo = sesData.taskInfo;
    
    %% Setup

    dfParams.kFoldIter = 10;
    dfParams.kFold = 5;
    dfParams.xValVerbosity = 0;
%    dfParams.selectFeaturesMinMax = [100 400];
%    dfParams.selectFeaturesMinMax = [10 1344]; %test using feat sel but no MRMR
    dfParams.ldaRegularization = 0.5; 
    
    
    % HARD CODED due to 2020.07.24 having some trials that did not match.
    % uTrialDuration = unique(taskInfo.trialDurationRounded);
    if isfield(sweepInfo,'uTrialDuration')
        uTrialDuration = sweepInfo.uTrialDuration;
    else
        uTrialDuration = [50 100, 200];
    end
    
    
    allLabels = {};
    allOffsets = {};
    
    
    isLprintfAvailable = ~isempty( which('lprintf') );
    if isLprintfAvailable
        % Because @fprintf below, matlab thinks this is now a variable,
        % so, we must set it.
        lprintf = @lprintf; 
    else
        lprintf = @fprintf;
    end
    

    attemptWindow = sweepInfo.attemptWindow;
    attemptOffset = sweepInfo.attemptOffset;
    noAttemptWindow = sweepInfo.noAttemptWindow;
    enforceTransitionMatrix = sweepInfo.enforceTransitionMatrix;

    selGest = ismember(taskInfo.labels,gestures); % Select trials for only certain gestures

    if ~strcmp(sweepInfo.labelType,'gestures')
        taskInfo.labels(:) = 'Attempt';
    end

    stepSize = sweepInfo.stepSize;
    featInds = sweepInfo.featInds;
    
    
    %% Sweep loop
    for di = 1:length(uTrialDuration)

        dur = uTrialDuration(di);

        offsets = attemptOffset(1):stepSize:(dur+attemptOffset(2));
        allOffsets{di} = offsets;
        clear allPerfOffsets

        % Select Trials
        selDur = ismember(taskInfo.trialDurationRounded, dur); % Select all trials for a given duration

        selTrl = selDur & selGest;

        events = taskInfo.startStops(selTrl, 1);

        if isempty(noAttemptWindow)
            noActionStStps = taskInfo.noActionStartStops;
        else
            noActionStStps = taskInfo.noActionStartStops(:,1) + noAttemptWindow;
        end
        
        noActionStStps(~selGest,:) = []; %remove no action periods from non-selected gesture trials (because T5 did not have extra no action trials)

        % It appears that what is being done here is removing the no action epochs from the events that are being used in the decoder
        % because the window over which we are evaluating the LDA (e.g. up to 2 second after end of action) overlaps with the period 
        % we are defining as our "no_action" period ([25 150], or 0.5-3sec after end of action).
        % The result is that for the LDA accuracy of the 1s duration attempts, we are using the no_action periods after the trials
        % of the 2s and 4s attempts.
        rmvNoAction = any(ismember(noActionStStps,RowColon(events + attemptOffset)),2);
        stopStStps = noActionStStps(~rmvNoAction,:);

        
        if isLprintfAvailable
            lprintf(); % Init print buffer / counter
        end
        
        for oi = 1:length(offsets)
            offset = offsets(oi);
            lprintf('%.0f sec attempts %03d offset   %02d of %02d (%.1f)\n', dur./50, offset, oi, length(offsets), oi/length(offsets)*100)



            % Select trials. Append no_action trials
            labels = cat(1, taskInfo.labels(selTrl), repmat({'no_action'},size(stopStStps,1),1));
            labels = removecats(labels);
            eventStStps = events(:) + attemptWindow + offset;


            startStops = cat(1, eventStStps, stopStStps);
            [~,si] = sort(startStops(:,1));

            labels = labels(si);
            startStops = startStops(si,:);

            %% Sanity check trial length
            trialLen = diff(startStops,[],2)+1;
            if any(trialLen<1)
                error('Invalid trial length')
            end
            rmvThresh = median(trialLen)*2;
            rmv = (trialLen > rmvThresh);
            if any(rmv)
                if isLprintfAvailable
                    lprintf('append-lprintf','Removed %d of %d trials (trial len > %d)\n',sum(rmv), length(rmv), rmvThresh)
                else
                    lprintf('Removed %d of %d trials (trial len > %d)\n',sum(rmv), length(rmv), rmvThresh)
                end
            end

            trialLen(rmv) = [];
            labels(rmv) = [];
            startStops(rmv,:) = [];
            
            
            %% Crossvalidate with LDA

            params     = MergeBintoA(dfParams, []);
            trainInds  = RowColon(startStops);
            flatLabels = GetLabelsForStartStops(labels,startStops);
            %         cu(flatLabels)
            
            % Create decoder instance
            trialBounds = cumsum( [1; trialLen] );
            params.trialBoundaries = trialBounds;
            m = dMultistate(params, feat(trainInds,featInds), flatLabels);
            m.params.smoothDataAmount = [5 0];
            m.params.enforceTransitionMatrix = enforceTransitionMatrix;
            m.xVal();
%             disp(sum(m.info.featureSelection.inds <193));disp(sum(m.info.featureSelection.inds <385));disp(sum(m.info.featureSelection.inds >385));

%             feat2 = smoothdata(feat,"movmean",5);
%%      
            useBuiltInLDA = 0;
            if useBuiltInLDA
                recall = GetRecallwithBuiltInLDA( feat(trainInds,featInds),flatLabels, params);
                m.performance.recall = recall;
            end
            
            allPerfOffsets(oi) = m.performance;

%%
        end
        allLabels{di} = labels;
        allPerfOffsetCell{di} = allPerfOffsets;
    end


    %% Output
    swp.perf = allPerfOffsetCell;
    swp.offsets = allOffsets;
    swp.durations = uTrialDuration;
    swp.info = sweepInfo;
    swp.info.labels = allLabels;
    swp.info.blocks = sesData.blocks;
    m.ClearDataAndLabels();
    swp.decoder = m;

end


function recall = GetRecallwithBuiltInLDA(X,Y,params)

% Define number of folds for cross-validation
numFolds = params.kFold;
iterations = round(params.kFoldIter / numFolds);

for i = 1:iterations
% Create a cross-validation partition
cv = cvpartition(Y, 'KFold', numFolds);

% Initialize recall storage for each fold
recallPerFold = zeros(numFolds, numel(unique(Y))); % Rows: folds, Columns: classes

% Perform cross-validation
for fold = 1:numFolds
    % Training and test indices for this fold
    trainIdx = training(cv, fold);
    testIdx = test(cv, fold);
    
    % Train the LDA model on training data
    ldaModel = fitcdiscr(X(trainIdx, :), Y(trainIdx));
    
    % Predict on the test data
    YTest = Y(testIdx);
    YPred = predict(ldaModel, X(testIdx, :));
    
    % Compute confusion matrix
    confMat = confusionmat(YTest, YPred, 'Order', unique(Y));
    
    % Calculate recall for each class
    % Recall = True Positives / (True Positives + False Negatives)
    recall = diag(confMat) ./ sum(confMat, 2);
    
    % Store recall for this fold
    recallPerFold(fold, :) = recall';
end

% Average recall across folds for each class
meanRecall(:,i) = mean(recallPerFold, 1);

% Display recall for each fold and mean recall
% disp('Recall for each fold:');
% disp(array2table(recallPerFold, 'VariableNames', unique(Y)));
% disp('Average recall across folds:');
% disp(array2table(meanRecall, 'VariableNames', unique(Y)));
end
recall = mean(meanRecall,2);
end