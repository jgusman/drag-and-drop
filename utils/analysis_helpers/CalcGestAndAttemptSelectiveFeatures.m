function [kwp, kwp_allDur, LDAacc, LDAacc_allDur] = CalcGestAndAttemptSelectiveFeatures(sesData,selectivityTypes,gestures,getBasicLDA)

if nargin < 3
    gestures = unique(sesData.taskInfo.labels);
end

if nargin < 4
    getBasicLDA = false;
end

%% Average features for each trial start stop
taskInfo = sesData.taskInfo;

features = sesData.feat;
nFeat = size(features,2);

durs = taskInfo.trialDurationRounded;
[cnt_uDurations,uDurations] = hist(durs,unique(durs));

uDurations(cnt_uDurations<3) = [];  %remove one off durations (or at least those with less than 3 trials (arbitrary) )

outlierThreshold = 5;

kwp = ones(nFeat,length(selectivityTypes),length(uDurations));
kwp_allDur = ones(nFeat,length(selectivityTypes));
LDAacc = zeros(round(nFeat/192),length(uDurations),length(selectivityTypes));
LDAacc_allDur = zeros(round(nFeat/192),length(selectivityTypes));


%% Get outlier trials
for jj = 1:length(uDurations)
    % Select all trials for a given duration
    selDur = ismember(taskInfo.trialDurationRounded, uDurations(jj));

    % Select trials only certain gestures
    selGest = ismember(taskInfo.labels,gestures);

    % Get average outlier per trial to exclude trials
    outlierWin = 1:uDurations(jj);
    outlierPerTrial = mean(GetEpochOfData2( taskInfo.prctNS5Outliers, taskInfo.startStops(:,1), outlierWin,'edt' ),3);
    excludeTrials = outlierPerTrial(:) > outlierThreshold;
    selTrlEachDur(:,jj) = selDur & selGest & ~excludeTrials;
    fprintf('Percent of %ds Trials Excluded: %0.0f%%\n', uDurations(jj)/50,mean(excludeTrials)*100)

    outlierWin_BL = 1:mode(diff(taskInfo.noActionStartStops,1,2)); %using whole no action period for no-action data
    outlierPerTrial_BL = mean(GetEpochOfData2( taskInfo.prctNS5Outliers, taskInfo.noActionStartStops(:,1), outlierWin_BL,'edt' ),3);
    excludeTrials_BL = outlierPerTrial_BL(:) > outlierThreshold;
    selTrlEachDur_BL(:,jj) = selGest & ~excludeTrials_BL; %using no action epochs from the same durations
    
    fprintf('Percent of Baseliine Trials Excluded: %0.0f%%\n', mean(excludeTrials_BL)*100)
end    
    

for sT = 1:length(selectivityTypes)
    selectivityType = selectivityTypes{sT};

    isAttemptSel = strcmp(selectivityType,'attempt');


    %% Calculate on ALL trials (i.e. all durations)
    selTrlAll = any(selTrlEachDur,2); % all trials with trial outliers excluded
    selTrlAll_BL = any(selTrlEachDur_BL,2); % all trials with trial outliers excluded

    % Get trial starts and stops and labels
    trStStp = taskInfo.startStops(selTrlAll,:);

    if isAttemptSel
        trStStp_BL = taskInfo.noActionStartStops(selTrlAll_BL,:);
        trStStp = cat(1,trStStp,trStStp_BL);
        labels = categorical( [ones(sum(selTrlAll),1); zeros(sum(selTrlAll_BL),1)], [0,1],{'Attempt','no_action'});
    else   %if gesture selective       
        labels = taskInfo.labels(selTrlAll);
    end

    [kwp_allDur(:,sT),LDAacc_allDur(:,sT)] = TestKW(trStStp,labels,features,getBasicLDA);


    %% Calculate on trial durations individually
    for jj = 1:length(uDurations)
    
        % Get trial starts and stops and labels
        selTrl = selTrlEachDur(:,jj);
        trStStp = taskInfo.startStops(selTrl,:);
    
        if isAttemptSel

            selTrl_BL = selTrlEachDur_BL(:,jj); %using no action epochs from the same durations
            trStStp_BL = taskInfo.noActionStartStops(selTrl_BL,:);
    
            trStStp = cat(1,trStStp,trStStp_BL);
            labels = categorical( [ones(sum(selTrl),1); zeros(sum(selTrl_BL),1)], [0,1],{'Attempt','no_action'});
        else   %if gesture selective       
            labels = taskInfo.labels(selTrl);
        end
    
        [kwp(:,sT,jj), LDAacc(:,sT,jj)] = TestKW(trStStp,labels,features,getBasicLDA);
        
    end

end

end


function [kwp,LDAacc] = TestKW(trStStp,labels,features,getBasicLDA)

        nTrials = size(trStStp,1);
        nFeat = size(features,2);
        avgData = nan(nTrials,nFeat);
        
        for ii = 1:nTrials
            inds = trStStp(ii,1):trStStp(ii,2);
            avgData(ii,:) = mean(features(inds,:));
        end
    
        % Compute Kruskal-Wallis p value per feature
        for f = 1:nFeat
            d = avgData(:,f); %average feature value per trial
            kwp(f) = kruskalwallis(d(:), labels(:), 'off');
        end
        
        % optionally get LDA accuracies using each set of feature types (e.g. LDA accuracy when just TX feats are used, or 20-39Hz LFPs are used, etc.)
        if getBasicLDA
            for fT = 1:round(nFeat/192)
                featInds = ((fT-1)*192+1) : (fT*192);
                LDA_model = fitcdiscr(avgData(:,featInds),labels);
                cvModel = crossval(LDA_model, 'KFold', 10); % 10-fold cross-validation
                cvLoss = kfoldLoss(cvModel); % Misclassification rate
                LDAacc(fT) = (1 - cvLoss) * 100;
        %         fprintf('Cross-Validated Accuracy: %.2f%%\n', LDA_Accuracy);
            end
        else
            LDAacc = zeros(round(nFeat/192),1);
        end

end