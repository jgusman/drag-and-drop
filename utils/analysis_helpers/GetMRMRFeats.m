
function [gestSelFeat,gestTopFeat,attemptSelFeat,attemptTopFeat] = GetMRMRFeats(sesData,kwp,numTopFeats,dur)

if nargin < 3
    numTopFeats = 400;
end
if nargin < 4
    dur = 200;
end

[nFeat,~,nDur] = size(kwp);

uDurations = [50 100 200];

selDur = uDurations == dur;



    trStStp = sesData.taskInfo.startStops;
    labels = sesData.taskInfo.labels;

    gestSelFeat = find(kwp(:,1,selDur) < 0.001); % gest;
    

        % MRMR on each training step
        trainInds = RowColon(trStStp);
        labelPerStep = [];
        for lbl = 1:length(labels)
            labelPerStep = [labelPerStep; repmat(labels(lbl),length(trStStp(lbl,1):trStStp(lbl,2)),  1 )];
        end
        [idx,scores] = fscmrmr(sesData.feat(trainInds,gestSelFeat),labelPerStep);
%         if numTopFeats > length(idx)
%             numTopFeats_gest = round((2/3)*idx);
%         else
            numTopFeats_gest = numTopFeats;
%         end
        gestTopFeat = gestSelFeat(idx(1:numTopFeats_gest)); %select top 400
    %     % MRMR on trial averages
    %     trialMeans = [];
    %     for ii = 1:length(labels)
    %         trialMeans(ii,:) = mean(sesData(sesI).feat(trStStp(ii,1):trStStp(ii,2),:));
    %     end
    %     [idx2,scores2] = fscmrmr(trialMeans(:,gestSelFeat),labels);
    % 



    attemptSelFeat = find(kwp(:,2,selDur) < 0.001); % attempt;


        % MRMR on each training step
        trainInds = RowColon(trStStp);
        noActionStStp = sesData.taskInfo.startStops;
    
        noActionTrainInds = RowColon(sesData.taskInfo.noActionStartStops);
    
        labelPerStep = ones(1,length(trainInds));
        labelPerStep = [labelPerStep, zeros(1,length(noActionTrainInds))];
    
        [idx,scores] = fscmrmr(sesData.feat([trainInds noActionTrainInds],attemptSelFeat),labelPerStep);
%         if numTopFeats > length(idx)
%             numTopFeats_attempt = round((2/3)*idx);
%         else
            numTopFeats_attempt = numTopFeats;
%         end
        attemptTopFeat = attemptSelFeat(idx(1:numTopFeats_attempt)); %select top 400


    



end