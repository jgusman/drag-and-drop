function [predState,latchState, gLiks, lLiks] = PredictFromLatchDecoder(dGesture,dLatch,features, opt)

if nargin == 0
    fprintf('\nTesting the predict logic...\n')
    TestLocalPredictLogic();
    fprintf('    Testing passed!\n')
    if nargout
        [predState,latchState,gLiks, lLiks] = deal([]);
    end
    return
end

if nargin < 4, opt = []; end
dfOpt.minGestureBeforeLatch_steps = 20;
dfOpt.gesturesThreshold = []; % if empty use mean( prctile(gLiks,95) ); %0.85;
dfOpt.latchThreshold = 0.9;
dfOpt.no_actionState = 1;
opt = MergeBintoA(dfOpt,opt);


%% Predict from features
if isa(dGesture, 'dMultistate')
    gLiks = dGesture.Predict(features);
    lLiks = dLatch.Predict(features);
    % Only care about the is-latch likelihood
    lLiks = lLiks(:,2);
elseif ismatrix(dGesture)
    removeGDims = sum(dGesture)==0;
    removeLDims = sum(dLatch)==0;
    
    gLiks = dGesture(:,~removeGDims);
    lLiks = dLatch(:,~removeLDims);
    
end

if isempty(opt.gesturesThreshold)
    ignoreOLInds = sum(gLiks,2)~=0;
    opt.gesturesThreshold = mean( prctile(gLiks(ignoreOLInds,:),95) );
    
elseif opt.gesturesThreshold > 1
    ignoreOLInds = sum(gLiks,2)~=0;
    opt.gesturesThreshold = mean( prctile(gLiks(ignoreOLInds,2:end),opt.gesturesThreshold) );
end
disp(opt.gesturesThreshold)

%% Loop through each step and apply thresholding/latching

[predState,latchState] = LocalPredictLogic(gLiks, lLiks, opt);

%%


end



function [stateDecode,latchState,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt)

nPoints = size(lLiks,1);

isLatched = false;
IsNoActionState = @(xx) xx == opt.no_actionState;
consecutiveDecodeCount = 0;
previousDecodeState = opt.no_actionState;

stateDecode = repmat(opt.no_actionState,nPoints,1);
latchState = false(nPoints,1);

for ii = 1:nPoints
    [largestGestLikVal,gestPredState] = max(gLiks(ii,:));
    isLatchAboveThresh = lLiks(ii) > opt.latchThreshold;
    isGestureAboveThresh = largestGestLikVal > opt.gesturesThreshold;
    
    if ~isGestureAboveThresh
        gestPredState = opt.no_actionState;
    end
    
    % This check should be before increment decodes if we only want to
    % increment when not latched...
    if isLatched
%         isLatchAboveThresh = ~(lLiks(ii) < 0.1);
        if ~isLatchAboveThresh
            isLatched = false;
        end
    end
    
    
    % Increment consecutive decodes
    if isGestureAboveThresh && gestPredState == previousDecodeState && ~IsNoActionState(gestPredState)
        % Do we want to increment currentDecodeState when latched?
        consecutiveDecodeCount = consecutiveDecodeCount + 1;
    else
        consecutiveDecodeCount = 0;
    end
    
    
    % Should be last
    if ~isLatched
        previousDecodeState = gestPredState;
        
        % Should we set the latch?
        if (consecutiveDecodeCount > opt.minGestureBeforeLatch_steps && ...
            ~IsNoActionState(previousDecodeState) && ...
            isLatchAboveThresh)
        
            isLatched = true;
        end
    end
    
    
    
    stateDecode(ii) = previousDecodeState;
    latchState(ii) = isLatched;
end

end



%% Testing the logic

function TestLocalPredictLogic()
    opt.minGestureBeforeLatch_steps = 5;
    opt.gesturesThreshold = 0.8;
    opt.latchThreshold = 0.2;
    opt.no_actionState = 1;


    %% Check that below/above threshold gesture liks output the right state
    nTime = 5;
    nStates = 3;
    [gLiks, lLiks] = InitLiks(nTime,nStates);
    gLiks(nTime-1:nTime,3) = 0.9;

    [stateDecode,isLatched,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt);

    isCorrect = IsCorrect(stateDecode,[1 1 1 3 3], ...
                    isLatched(end),0, ...
                    consecutiveDecodeCount, 1);

    if ~isCorrect
        error('Gesture likelihood decode check failed')
    end


    %% Check if latch happens
    [gLiks, lLiks] = InitLiks(5+opt.minGestureBeforeLatch_steps,nStates);
    gLiks(end-(opt.minGestureBeforeLatch_steps+1):end,nStates) = 0.9;
    lLiks(:) = 0.9;
    [stateDecode,isLatched,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt);

    expectedState = ones(size(lLiks));
    expectedState(4:end) = nStates;
    [isCorrect, isCorState,isCorLatch,isCorConsec] = IsCorrect(stateDecode,expectedState, ...
            isLatched(end),1, ...
            consecutiveDecodeCount, opt.minGestureBeforeLatch_steps+1);
    if ~isCorrect
        error('Latch enable check failed')
    end

    %% Check that latch works

    [gLiks, lLiks] = InitLiks(5+opt.minGestureBeforeLatch_steps+5,nStates);
    gLiks(4:4+opt.minGestureBeforeLatch_steps+1,nStates) = 0.9;
    lLiks(:) = 0.9;
    [stateDecode,isLatched,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt);

    expectedState = ones(size(lLiks));
    expectedState(4:end) = nStates;
    [isCorrect, isCorState,isCorLatch,isCorConsec] = IsCorrect(stateDecode,expectedState, ...
            isLatched(end),1, ...
            consecutiveDecodeCount, 0);
    if ~isCorrect
        error('Latch check failed')
    end

    %% Check that latch resets

    [gLiks, lLiks] = InitLiks(5+opt.minGestureBeforeLatch_steps+5,nStates);
    gLiks(4:4+opt.minGestureBeforeLatch_steps+1,nStates) = 0.9;
    lLiks(:) = 0.9;
    lLiks(end-1:end) = 0;

    [stateDecode,isLatched,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt);

    expectedState = ones(size(lLiks));
    expectedState(4:end-2) = nStates;
    [isCorrect, isCorState,isCorLatch,isCorConsec] = IsCorrect(stateDecode,expectedState, ...
            isLatched(end),0, ...
            consecutiveDecodeCount, 0);
    if ~isCorrect
        error('Latch reset check failed')
    end


    %% Check a decoded gesture does not override the latch

    [gLiks, lLiks] = InitLiks(5+opt.minGestureBeforeLatch_steps+5,nStates);
    gLiks(4:4+opt.minGestureBeforeLatch_steps+1,nStates) = 0.9;
    lastValidDecIdx = 4+opt.minGestureBeforeLatch_steps+1;
    gLiks(lastValidDecIdx+1,2) = 0.9;
    lLiks(:) = 0.9;
    lLiks(end-1:end) = 0;

    [stateDecode,isLatched,consecutiveDecodeCount] = LocalPredictLogic(gLiks, lLiks, opt);

    expectedState = ones(size(lLiks));
    expectedState(4:end-2) = nStates;
    [isCorrect, isCorState,isCorLatch,isCorConsec] = IsCorrect(stateDecode,expectedState, ...
            isLatched(end),0, ...
            consecutiveDecodeCount, 0);
    if ~isCorrect
        error('Latch-gesture override check failed')
    end




end

function [gLiks, lLiks] = InitLiks(nTime,nStates)
gLiks = zeros(nTime,nStates);
lLiks = zeros(nTime,1);
end

function [isCorrect, isCorState,isCorLatch,isCorConsec] = IsCorrect(stateDecode,expectedState, isLatched, expectedLatch, consecutiveDecodeCount, expectedConsecutiveDecodeCount)
isCorState = isequal(stateDecode(:),expectedState(:));
isCorLatch = isequal(isLatched(:),expectedLatch(:));
isCorConsec = isequal(consecutiveDecodeCount(:),expectedConsecutiveDecodeCount(:));
isCorrect = isCorState && ...
            isCorLatch && ...
            isCorConsec;
end
