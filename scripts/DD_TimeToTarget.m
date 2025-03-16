%% Time to Target Plots (Figs 6C-D)
% T11's Closed Loop performance of Drag and Drop task as the time to target duration
% for the "Center Out" and "Return" stages of Drag trials for Session 1 (6C) and Session 2 (6D)

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------


% set figure save folder
ddFigHelper.ResetParams()
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'DD_TimeToTarget'))

%% Calculate Trial Durations
sPerformance = GetPerformance(sesData_DD_T11);

%% Plot Figures 6(C) and 6(D): Drag and Drop Trial Stage Durations / Time to target
metric = 'moveTime_sec'; 
fH = PlotPerformanceBar(sPerformance,metric);

% Save Plots
legend(fH(1).CurrentAxes,'off') %turn off legend
ddFigHelper.SaveFigure(fH(1),'Figure 6C - DD Time2Targ Comparison - Session 1')

legend(fH(2).CurrentAxes,'off') %turn off legend
ddFigHelper.SaveFigure(fH(2),'Figure 6D - DD Time2Targ Comparison - Session 1')

%% Log Statistics
% Below creates a log file in the saveDir/logSubDir/Trial Duration Statistics- Move vs Click vs Drag.txt
ddFigHelper.LogPrints('Time2Targ Statistics - Move vs Click vs Drag');
PrintPerformanceStatistics(sPerformance,metric)
ddFigHelper.LogPrints(); % This turns off the log






%% %%%%% HELPER FUNCTIONS %%%%% %%


function sPerformance = GetPerformance(sesData)
nSessions = length(sesData);
fprintf('\n\nTrial duration analysis\n');

for sesI = 1:nSessions
    taskInfo = sesData(sesI).taskInfo;
    
    [sPerf,analyzeBlocks,excludedTrialNums] = SelectTrialStages(taskInfo);
    sPerf = GetMoveTime(sPerf,taskInfo);
    
    fprintf('Drag and Drop Session %s: Blocks %s\n', sesData(sesI).sessionNum, BlocksToString( analyzeBlocks ) )
    fprintf(' Excluding %d trials due to ns5 outliers\n\n',length(excludedTrialNums))
    sPerformance(sesI,:) = sPerf;
end

end



function [sPerf,analyzeBlocks,excludedTrialNums] = SelectTrialStages(taskInfo)

n = 1;

sPerf(n).trialType = 'Move Only';
sPerf(n).stageName = 'Center Out';
sPerf(n).stageStart = 'Prepare - Direction';
sPerf(n).epochOffsets = [1 2]; % (move, wait) From prepare epoch
n = n + 1;

sPerf(n).trialType = 'Move Only';
sPerf(n).stageName = 'Return';
sPerf(n).stageStart = 'Prepare - Direction';
sPerf(n).epochOffsets = [3 4]; % From prepare epoch
n = n + 1;

sPerf(n).trialType = 'Click';
sPerf(n).stageName = 'Center Out';
sPerf(n).stageStart = 'Prepare - Direction + Gesture click';
sPerf(n).epochOffsets = [1 2];
n = n + 1;

sPerf(n).trialType = 'Click';
sPerf(n).stageName = 'Return';
sPerf(n).stageStart = 'Prepare - Direction + Gesture click';
sPerf(n).epochOffsets = [3 4];
n = n + 1;

sPerf(n).trialType = 'Drag';
sPerf(n).stageName = 'Center Out';
sPerf(n).stageStart = 'Prepare - Direction + Gesture'; %i.e. first stage of trial
sPerf(n).epochOffsets = [1 2];
n = n + 1;

sPerf(n).trialType = 'Drag';
sPerf(n).stageName = 'Return';
sPerf(n).stageStart = 'Prepare - Direction + Gesture';
sPerf(n).epochOffsets = [3 4];

%% Exclude Trials
% exclude data from calibration trials
% exclude data with excessive NS5 outliers

isExcludedCalTrial =  ~taskInfo.isKinDecoderCL | taskInfo.kinErrorAttenuation ~= 0;
isExcludedOutliersTrial = ddHelper.ExcludePerformanceTrialsWithNS5Outliers(taskInfo);
isExcludedTrial = isExcludedCalTrial | isExcludedOutliersTrial;
blockNum = taskInfo.blockNumber;

for mtI = 1:length(sPerf)
    isSelEpochsPreExlude = (ismember(taskInfo.trialStage,sPerf(mtI).stageStart));
    isSelEpochs = isSelEpochsPreExlude & ~isExcludedTrial;
    
    selEpochs = find(isSelEpochs);
    epochStartStops = selEpochs + sPerf(mtI).epochOffsets;

    sPerf(mtI).epochStartStops = epochStartStops;
    sPerf(mtI).isSelEpochsPreExlude = isSelEpochsPreExlude;
    sPerf(mtI).blockNumber = taskInfo.blockNumber(isSelEpochs);
end

% For info
stageEpochLogical = false(size(isExcludedOutliersTrial));
for ii = 1:size(sPerf)
    stageEpochInds = find( sPerf(ii).isSelEpochsPreExlude & ~isExcludedCalTrial ) + sPerf(ii).epochOffsets;
    stageEpochInds = RowColon(stageEpochInds);
    stageEpochLogical(stageEpochInds) = true;
end

excludedOutlierEpochs = stageEpochLogical & isExcludedOutliersTrial;
excludedTrialNums = unique(taskInfo.trialNum(excludedOutlierEpochs));
analyzeBlocks = unique(taskInfo.blockNumber(~isExcludedTrial));

end


function sPerf = GetMoveTime(sPerf,taskInfo)

% Find move times
for mtI = 1:length(sPerf)

    epochStartStops = sPerf(mtI).epochStartStops;
    durationInds = nan(size(epochStartStops));

    % Stop at target contact / last target contact
    durationInds(:,1) = taskInfo.startStops(epochStartStops(:,1),1);
    for ii = 1:size(epochStartStops,1)
        onTargStStps = taskInfo.onTargetStartStops{epochStartStops(ii,1)};
        if isempty(onTargStStps)
            % use the next stage start
            onTargStStps = taskInfo.startStops(epochStartStops(ii,2),1);
        else
            if onTargStStps(end,1) > taskInfo.startStops(epochStartStops(ii,2),1)  %is end of ontarget period later than beginning of Wait stage? (debug i think)
                fprintf('Large onset\n')
            end
        end
        %                 durationInds(ii,2) = onTargStStps(1,1); % First on target
        durationInds(ii,2) = onTargStStps(end,1); % last on target (i.e. first time step of the last time target was reached
        % i.e. we dont count time target was traversed over quickly)
    end
        
    moveSteps = diff(durationInds,[],2);
    moveTime = moveSteps./50;

    sPerf(mtI).durationStStp = durationInds;
    sPerf(mtI).moveTime_sec = moveTime;
end

end




function [fh,sPerformance] = PlotPerformanceBar(sPerformance,metric)

nSessions = size(sPerformance,1);

[trialTypes,~,trialTypeInds] = unique({sPerformance(1,:).trialType},'stable'); %  {'Drag','Click','MoveOnly'}; %should match the num of columns in compareInds
[stageNames,~,stageInds] = unique({sPerformance(1,:).stageName},'stable'); %{'Center Out', 'Return'}; %should match the num of rows in compareInds

nTrialTypes = length(trialTypes);
nStages = length(stageNames);
nGroups = size(sPerformance,2); % each group is set of trials for each trial type and stage type (e.g. Return/MoveOnly)

data = nan([nSessions,nGroups,100]);  %100 arbitrary large (greater than most number of trials for a group)

trialTypeColors = ddHelper.trialTypeColors;

% Collect Data
for sesI = 1:nSessions
    for i = 1:size(sPerformance,2)
        d = sPerformance(sesI,i).(metric);
        data(sesI,i,1:length(d)) = d;     
    end
end

% Plot each session on its own
for sesI = 1:nSessions
    % Get Medians
    dataMedians = nanmedian(data(sesI,:,:),3);
    % Get CIs
    dataCIs = bootci(1000,@nanmedian,squeeze(data(sesI,:,:))');  % 95% CI of the median
    dataCIs = dataCIs-dataMedians;
    
    % Reorganize for plotting
    fh(sesI) = ddFigHelper.CreateFigure(400+sesI);
    clf
    
    interPlotDist = 0.3;
    interBarDist = 0.05;

    barWidth = (1-interPlotDist-interBarDist*(nTrialTypes-1))  / nTrialTypes;
    oi = (1-interPlotDist-barWidth) / 2; %offset

    xs = linspace(-oi,oi,nTrialTypes);

    xtickVals = 1:nStages;

    for ii = 1:nTrialTypes
        cmpI = trialTypeInds == ii;
        
        % Plot Bar Plots
        bX = xtickVals + xs(ii);
        bh(ii) = bar(bX,dataMedians(cmpI));
        bh(ii).BarWidth = barWidth;
        bh(ii).FaceColor = trialTypeColors(ii,:); 
        
        % Plot Error Bars
        errX = bh(ii).XEndPoints;% ii + xs(jj); %same as bx?
        eb = errorbar(errX,dataMedians(cmpI),dataCIs(1,cmpI),dataCIs(2,cmpI),'.k', 'MarkerEdgeColor', 'k', 'Color','k');
        eb.CapSize = 2; %width of top and bottom bar things

    end
    legend(bh,trialTypes)

    ax = gca;
    ax.XTick = xtickVals;
    ax.XTickLabel = stageNames;
    ax.FontSize = 8;
    switch metric
        case 'moveTime_sec'
            ylabel('Time to Target (s)','FontSize',7)
            xlim([min(xtickVals)-(interPlotDist+barWidth), max(xtickVals)+(interPlotDist+barWidth)])
            ylim([0 6])
        case 'angleError'
            ylabel('Angular Error (°)')
        otherwise
    end

    ddFigHelper.UpdateFigDims([2.5 1])
end

end



function PrintPerformanceStatistics(sPerformance,metric)
    
    nSessions = size(sPerformance,1);
    
    [trialTypes,~,trialTypeInds] = unique({sPerformance(1,:).trialType},'stable'); %  {'Drag','Click','MoveOnly'}; %should match the num of columns in compareInds
    [stageNames,~,stageInds] = unique({sPerformance(1,:).stageName},'stable'); %{'Center Out', 'Return'}; %should match the num of rows in compareInds
    
    nTrialTypes = length(trialTypes);
    nStages = length(stageNames);
    nGroups = size(sPerformance,2); % each group is set of trials for each trial type and stage type (e.g. Return/MoveOnly)
    
    data = nan([nSessions,nGroups,100]);  %100 arbitrary large (greater than most number of trials for a group)
        
    % Collect Data
    for sesI = 1:nSessions
        for i = 1:size(sPerformance,2)
            d = sPerformance(sesI,i).(metric);
            data(sesI,i,1:length(d)) = d;     
        end
    end
    
    % Statistics for each session
    for sesI = 1:nSessions
        
        for i = 1:nStages
            cInds = stageInds==i;
            d = squeeze(data(sesI,cInds,:))';
            
            alpha = 0.05;
            fprintf('\n\n~~~Session %d Statistics (%s)~~~\n',sesI,stageNames{i});
    
            for j = 1:nTrialTypes
    
                trialDuration = d(:,j);
                trialDuration(isnan(trialDuration)) = [];
                timeout_thresh = 25;
                isTimeout = (trialDuration-timeout_thresh) > 0;
    
                fprintf('\n%s\n',trialTypes{j})
                prtTimeout = sum(isTimeout)/length(isTimeout)*100;
                fprintf(' - Mean   %2.2f (std %.2f) seconds\n', mean(trialDuration), std(trialDuration))
                fprintf(' - Median %2.2f (iqr %.2f) seconds\n', median(trialDuration), iqr(trialDuration))
                fprintf(' - Percent timeout: %.2f%% (%d/%d)\n', prtTimeout, sum(isTimeout),length(isTimeout))  
            end
    
            KWandRSTests(d,trialTypes,alpha);
    
        end
    end

    % Statistics Both Sessions Combined 
    for i = 1:nStages
        cInds = stageInds==i;
        d1 = squeeze(data(1,cInds,:))'; %should functionalize this
        d2 = squeeze(data(2,cInds,:))';
        d = vertcat(d1,d2);

        alpha = 0.05;
        fprintf('\n\n~~~Sessions Combined Statistics (%s)~~~\n',stageNames{i});

        for j = 1:nTrialTypes

            trialDuration = d(:,j);
            trialDuration(isnan(trialDuration)) = [];
            timeout_thresh = 25;
            isTimeout = (trialDuration-timeout_thresh) > 0;

            fprintf('\n%s\n',trialTypes{j})
            prtTimeout = sum(isTimeout)/length(isTimeout)*100;
            fprintf(' - Mean   %2.2f (std %.2f) seconds\n', mean(trialDuration), std(trialDuration))
            fprintf(' - Median %2.2f (iqr %.2f) seconds\n', median(trialDuration), iqr(trialDuration))
            fprintf(' - Percent timeout: %.2f%% (%d/%d)\n', prtTimeout, sum(isTimeout),length(isTimeout))  
        end

        KWandRSTests(d,trialTypes,alpha);
    end
end