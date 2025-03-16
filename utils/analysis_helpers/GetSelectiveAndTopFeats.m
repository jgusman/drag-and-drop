
function [selFeat,topFeat] = GetSelectiveAndTopFeats(sesData,kwp,sigAlpha,numTopFeats,gestures)

if nargin < 3
    sigAlpha = 0.001;
end

if nargin < 4
    numTopFeats = 400;
end

if nargin < 5
    gestures = unique(sesData.taskInfo.labels);
end

uTrialDurations = [50 100 200];

[nFeat,~,nDur] = size(kwp);


selFeat = false(size(kwp));
topFeat = false(size(kwp));

for di = 1:nDur

    selGest = ismember(sesData.taskInfo.labels,gestures);
    selDur = sesData.taskInfo.trialDurationRounded == uTrialDurations(di);
    selTrl = selGest & selDur;

    trStStp = sesData.taskInfo.startStops(selTrl,:);
    labels = sesData.taskInfo.labels(selTrl,:);

    % Gesture Selective Features
    gestSelFeat = find(kwp(:,1,di) < sigAlpha); % gest;
    
    % MRMR on each training step
    trainInds = RowColon(trStStp);
    labelPerStep = [];
    for lbl = 1:length(labels)
        labelPerStep = [labelPerStep; repmat(labels(lbl),length(trStStp(lbl,1):trStStp(lbl,2)),  1 )];
    end
    [idx,scores] = fscmrmr(sesData.feat(trainInds,gestSelFeat),labelPerStep);

    numTopFeats_gest = numTopFeats;
    try
        gestTopFeat = gestSelFeat(idx(1:numTopFeats_gest)); %select top 400
    catch
        gestTopFeat = gestSelFeat; % if already fewer than 400 features in analysis (e.g. if just the 384 TX and SP feats)
    end

    selFeat(gestSelFeat,1,di) = true;
    topFeat(gestTopFeat,1,di) = true;


    % Attempt Selective Features
    attemptSelFeat = find(kwp(:,2,di) < sigAlpha); % attempt;

    % MRMR on each training step
    noActionStStp = sesData.taskInfo.noActionStartStops(selTrl,:);

    noActionTrainInds = RowColon(noActionStStp);

    labelPerStep = ones(1,length(trainInds));
    labelPerStep = [labelPerStep, zeros(1,length(noActionTrainInds))];

    [idx,scores] = fscmrmr(sesData.feat([trainInds noActionTrainInds],attemptSelFeat),labelPerStep);

    numTopFeats_attempt = numTopFeats;

    try
        attemptTopFeat = attemptSelFeat(idx(1:numTopFeats_attempt)); %select top 400
    catch
        attemptTopFeat = attemptSelFeat;
    end

    selFeat(attemptSelFeat,2,di) = true;
    topFeat(attemptTopFeat,2,di) = true;
    
end


end