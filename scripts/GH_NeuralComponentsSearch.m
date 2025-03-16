%% Neural Components Search (Figs 3C-D, S6, S7, and S8)
% Sweeps training multistate decoder with varying training window to look for onset and offset responses.
% Then plots the sweep results (using adjusted mcc).

% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------


% set figure save subfolder
ddFigHelper.ResetParams()
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'NeuralComponents'))

%% Define grid search resolution
increment = 5;  % grid search resolution; number of 20ms timestamps to increment 
                % in sweep of onset & offset times; default is 5 (100ms)

%% Figure 3(C): Grid search on T11 4s power grasp trials
[sweepResults119,fH119] = NCSearch(sesData_GH_T11,200,{'Right hand - power grasp'},[],increment,'Power, Both Sess');

%% Figure 3(D): Grid search on T5 4s power grasp trials
[sweepResults53,fH53] = NCSearch(sesData_GH_T5(1),200,{'Right hand - power grasp'},[],increment,'Power');

%% Figure S6: Grid search results (T11 and T5) for 1s, 2s, and 4s power grasp trials
[sweepResults117,fH117] = NCSearch(sesData_GH_T11,50,{'Right hand - power grasp'},[],increment,'Power, Both Sess');
[sweepResults118,fH118] = NCSearch(sesData_GH_T11,100,{'Right hand - power grasp'},[],increment,'Power, Both Sess');
% [sweepResults119,fH119] = NCSearch(sesData_GH_T11,200,{'Right hand - power grasp'},[],[],'Power, Both Sess'); % (uncomment if not previously calculated)

[sweepResults51,fH51] = NCSearch(sesData_GH_T5(1),50,{'Right hand - power grasp'},[],increment,'Power');
[sweepResults52,fH52] = NCSearch(sesData_GH_T5(1),100,{'Right hand - power grasp'},[],increment,'Power');
% [sweepResults53,fH53] = NCSearch(sesData_GH_T5(1),200,{'Right hand - power grasp'},[],[],'Power'); % (uncomment if not previously calculated)

%% Figure S7(A): Grid search results on T11 all (7-gestures) 4s trials - attempt vs no attempt decoding
[sweepResults1110,fH1110] = NCSearch(sesData_GH_T11,200,[],[],increment,'All Gest, Both Sess'); %binary decode (takes long time)

%% Figure S7(B): Grid search results on T5 all (3-gestures) 4s trials - attempt vs no attempt decoding
[sweepResults54,fH54] = NCSearch(sesData_GH_T5(1),200,[],[],increment,'All Gest'); %binary decode

%% Figure S8(A): Grid search results on T11 all (7-gestures) 4s trials - multi-gesture decoding
[sweepResults1111,fH1111] = NCSearch(sesData_GH_T11,200,[],true,increment,'All Gest, Both Sess'); % multistate decode (takes long time)

%% Figure S8(B): Grid search results on T5 all (3-gestures) 4s trials - multi-gesture decoding

[sweepResults55,fH55] = NCSearch(sesData_GH_T5(1),200,[],true,increment,'All Gest'); % multistate decode





%% %%% HELPER FUNCTIONS %%% %%

function [sweepResults,fH1,fH2] = NCSearch(sesData,duration,gestures,decodeMultistate,increment,appendStr)
    %
    clear swOpt
    
    swOpt.duration = duration;
    
    if ~isempty(gestures)
        swOpt.gestures = gestures;
    end
    if ~isempty(decodeMultistate)
        swOpt.decodeMultistate = decodeMultistate;
    end
    if ~isempty(increment)
        swOpt.increment = increment;
    end

    % Sweep Decoder
    sweepResults = TrainWinSweep(sesData,swOpt);
    
    % Plot MCC Sweep
    try    %using try statement so that if plotting breaks, will still output the sweepResults
        %%
        fH1 = PlotMCC2(sweepResults);
        
        % Save
        if sweepResults(1).swOpt.decodeMultistate
            decMultStr = 'Multistate';
        else
            decMultStr = 'Binary';
        end
        fileName = sprintf('%s Sweep, Dur %d, %s, %s',upper(sesData(1).participant),sweepResults(1).swOpt.duration,decMultStr,appendStr);
        ddFigHelper.SaveFigure(fileName);
        %%
        fH2 = PlotOnsetOffsetSustained(sweepResults);

        fileName = sprintf('%s Sweep, Dur %d, %s, %s COMPONENTS',upper(sesData(1).participant),sweepResults(1).swOpt.duration,decMultStr,appendStr);
        ddFigHelper.SaveFigure(fileName);
        
    catch ME
        fH1 = []; fH2 = [];
        warning('Plot or Save Failure:')
        warning(ME.message)
    end

end



function fH = PlotOnsetOffsetSustained(sweepResults)

    trainOpt = sweepResults(1).trainOpt;
    swOpt = sweepResults.swOpt;

    if length(sweepResults) > 1
        mcc = mean( cat(3,sweepResults(:).mccAdj), 3);
    else
        mcc = sweepResults.mccAdj;
    end

    mcc_on = mcc; %onset component
    mcc_off = mcc; %offset component
    mcc_sus = mcc; %sustained component
    mcc_on(:,trainOpt.trainStops>50) = 0;
    mcc_off(trainOpt.trainStarts<swOpt.duration-50,:) = 0;
    mcc_sus(:,trainOpt.trainStops<50) = 0;
    mcc_sus(trainOpt.trainStarts>swOpt.duration-50,:) = 0;

    onColor =  [     0    0.4470    0.7410];
    offColor = [0.8500    0.4       0.0980];
    susColor = [0.4940    0.1840    0.5560];

    fH = figure;
    tiledlayout(3,6)
    
    % Plot MCC adjusted
    nexttile([3 2])
    PlotMCC(gca,sweepResults,'adjusted');
    title('MCC_a_d_j')

    % Plot example trial likelihoods
    nexttile([1 4])
    [bestStart_off,bestStop_off] = PlotBestTrainWin(mcc_off,trainOpt,swOpt,offColor);
    title("offset",'Color',offColor)
    xlabel ''
    ylabel ''

    nexttile([1 4])
    [bestStart_sus,bestStop_sus] = PlotBestTrainWin(mcc_sus,trainOpt,swOpt,susColor);
    title("sustained",'Color',susColor)
    xlabel ''

    nexttile([1 4])
    [bestStart_on,bestStop_on] = PlotBestTrainWin(mcc_on,trainOpt,swOpt,onColor);
    title("onset",'Color',onColor)
    ylabel ''

    fH.Position(3) = 1350;
    
    % Plot squares on grid denoting local maxima
    nexttile(1)
    hold on
    plot(bestStart_on/50,bestStop_on/50,'s',"Color", onColor,"LineWidth",2.5,"MarkerSize",12,"ZData",1)
    plot(bestStart_off/50,bestStop_off/50,'s',"Color", offColor,"LineWidth",2.5,"MarkerSize",12,"ZData",1)
    plot(bestStart_sus/50,bestStop_sus/50,'s',"Color", susColor,"LineWidth",2.5,"MarkerSize",12,"ZData",1)

    AddColorBar()

end



function sweepResults = TrainWinSweep(sesData,swOpt)
% Sweeps eventStartWin and eventWinLen for attempt onset and offset
% 
% Output
%   sweepResults [nSes, 1]
%   sweepResults(sesI) = StructFromVars(allPerf, mccPerf, mccAdj, winSize, trainOpt);

% default options
dfOpt.featInds = 1:384;
dfOpt.sweepRangeStStp = [-75 75]; % how many steps before trial start and after trial stop to test
dfOpt.trialBufferStStp = [-150 100]; % how many steps before trial start and after trial stop to include as data in decoder
dfOpt.increment = 5;
dfOpt.decodeMultistate = false;  % multistate decoder or binary decoder where all gestures are treated as one "attempt" class
dfOpt.duration = 50;
dfOpt.gestures = [];

swOpt = MergeBintoA(dfOpt,swOpt);

trainStarts = swOpt.sweepRangeStStp(1) : swOpt.increment : swOpt.sweepRangeStStp(2)+swOpt.duration-swOpt.increment; %relative to event start
trainStops = swOpt.sweepRangeStStp(1)+swOpt.increment : swOpt.increment : swOpt.sweepRangeStStp(2)+swOpt.duration; %relative to event start

    tic
    for sesI = 1:length(sesData)
        taskInfo = sesData(sesI).taskInfo;
        feat = sesData(sesI).feat;

        fullTrialStStps = taskInfo.startStops + swOpt.trialBufferStStp;
        avgOutliers = [];
        for ii = 1:size(fullTrialStStps,1)
            avgOutliers(ii,1) = mean(taskInfo.prctNS5Outliers(RowColon(fullTrialStStps(ii,:))));
        end
        
        if isempty(swOpt.gestures)
            selTrials = taskInfo.trialDurationRounded == swOpt.duration; % using all gestures
        else
            selTrials = taskInfo.trialDurationRounded == swOpt.duration & ismember(taskInfo.labels, swOpt.gestures); 
        end

        selI = (avgOutliers(:) < 3) & selTrials;
        fullTrialStStps = fullTrialStStps(selI,:);

        trainInds = RowColon(fullTrialStStps);
        trialBounds = GetTrialBounds(fullTrialStStps);
        trainFeat = feat(trainInds,swOpt.featInds);

        trainLabels = GetLabelsForStartStops(taskInfo.labels(selI), fullTrialStStps);

        fprintf('\nSession %s (%d of %d)\n', sesData(sesI).date, sesI, length(sesData));

        % Sweep onset / offset
        clear allPerf
        mccPerf = nan(length(trainStarts),length(trainStops));
        winSize = nan(length(trainStarts),length(trainStops));
        
        lprintf();
        currentIter = 0;
        totalIters = length(trainStarts)*(length(trainStarts)+1)/2;

        startStops = taskInfo.startStops(selI,:);
        for startI = 1:length(trainStarts)
            for stopI = startI:length(trainStops)

                currentIter = currentIter+1;
                percentComplete = currentIter/totalIters*100;
                lprintf('(%0.2f/100) Starts: %02d of %02d | Stops: %02d of %02d | Estimated Time Remaining: %0.2f min.\n', percentComplete,startI,length(trainStarts), stopI - startI, length(startI:length(trainStops)), (toc * (1-percentComplete/100)/(percentComplete/100)) /60 );
                m = SetupDecoder(trainFeat,trainInds,trialBounds, startStops,trainLabels, trainStarts(startI), trainStops(stopI),swOpt.decodeMultistate);
                m.xVal();
                allPerf(startI,stopI) = m.performance;
                mccPerf(startI,stopI) = m.performance.mcc;
                winSize(startI,stopI) = trainStops(stopI)-trainStarts(startI);
            end
        end

        % Calculated Adjusted MCC
        mccAdj = mccPerf;
        uWinSize = unique(winSize);
        for w = 1:length(uWinSize)
            if uWinSize(w) == 0
                continue
            end
            winSizeInds = winSize==uWinSize(w);
            K = min(mccPerf(winSizeInds)); % "K is the minimum MCC achieved across all classifications with the same class imbalance(i.e.window size)."
            mccAdj(winSizeInds) = (mccPerf(winSizeInds) - K) / (1-K); % adjusted MCC formula as used in Dekleva et al, 2021
        end

        trainOpt = StructFromVars(trainFeat, trainInds, trialBounds, startStops, trainLabels, trainStarts, trainStops);

        sweepResults(sesI) = StructFromVars(allPerf, mccPerf, mccAdj, winSize, trainOpt,swOpt);

    end
    fprintf('\n\nDone! Elapsed Time: %.1f min\n',toc/60)
%     sound([sin(0:0.5:500),sin(0:0.6:1000)])

end


function m = SetupDecoder(feat,trainInds,trialBounds, startStops, labels, trainWinStart, trainWinStop, decodeMultistate)
%%
    % Get train inds 
%     trainWindows = startStops + [trainWinStart trainWinStop];
    trainWindows = startStops(:,1) + [trainWinStart trainWinStop]; %start and stop inputs are relative to event start

    taskTrainInds = RowColon(trainWindows);    
    % Map the specified train window to trainInds
    % Set the train window to be true
    % All other points are false
    trueTrain = ismember(trainInds,taskTrainInds(:))';

    if decodeMultistate %use Gesture labels
        % Gestures
        trainLabels = labels;
        trainLabels(~trueTrain) = 'no_action';
        trainLabels = removecats(trainLabels); %remove unused categories (in order to not confuse dMultistate)
        labelPerTrial = labels(trialBounds(1:end-1));
    else
        % Attempt vs no attempt
        trainLabels = categorical(double(trueTrain),0:1,{'no_action','attempt'});
        labelPerTrial = []; % Doesn't matter
    end

    m = dMultistate();
    m.params.trialBoundaries = trialBounds;
    m.params.labelPerTrial = labelPerTrial;
    m.data = feat;
    m.labels = trainLabels;
    
    m.params.enforceTransitionMatrix = 0; %for not using HMM, just LDA
%     m.params.smoothDataAmount = [5 0]; % mild smoothing of data (optional)
%     m.params.smoothLikAmount = [15 0]; % mild smoothing of likelihoods (optional)

end


function PlotMCC(ax,sweepResults,varargin)
% Plot MCC heatmap

    sampleRate = 50;
    swOpt = sweepResults(1).swOpt;

    useAdjustedMcc = false;
    if ~isempty(varargin)
        if strcmp(varargin{1},'adjusted')
            useAdjustedMcc = true;
        end
    end
    
    if length(sweepResults) == 1
        mccAdj = sweepResults.mccAdj;
        mccPerf = sweepResults.mccPerf;
    else
        mccAdj = mean( cat(3,sweepResults(:).mccAdj), 3);
        mccPerf = mean( cat(3,sweepResults(:).mccPerf), 3);
    end
    
    if useAdjustedMcc
        plotMcc = mccAdj';
    else
        plotMcc = mccPerf';
    end
    
    xLims = [swOpt.sweepRangeStStp(1), swOpt.duration+swOpt.sweepRangeStStp(2)] / sampleRate;
    yLims = [swOpt.sweepRangeStStp(1), swOpt.duration+swOpt.sweepRangeStStp(2)] / sampleRate;
    
    im = imagesc(ax,xLims, yLims, plotMcc);

    alphaData = ~isnan(plotMcc);
    im.AlphaData = alphaData; %remove non values
    
    ax.YDir = 'normal';
    ax.FontSize = 13;
    ax.LineWidth = 2;
    ax.FontWeight = 'bold';
    box off

    % set ticks and labels
    tickVals = -1:swOpt.duration/sampleRate+1;
    xticks(tickVals)
    yticks(tickVals)
    ax.XTickLabels([2 length(tickVals)-1]) = {'Att.', 'Rel.'};
    ax.YTickLabels([2 length(tickVals)-1]) = {'Att.', 'Rel.'};
    
    ax.XAxis.MinorTickValues = xLims(1):0.25:xLims(2);
    ax.XMinorTick = 'on';
    ax.YAxis.MinorTickValues = xLims(1):0.25:xLims(2);
    ax.YMinorTick = 'on';

    xlabel('Train Window Start (s)','FontSize',16)
    ylabel('Train Window Stop (s)','FontSize',16)

    ax.XLim = xLims;
    ax.YLim = yLims;
end


function fH = PlotMCC2(sweepResults)
% Plot MCC regular and MCC adjusted side-by-side
    fH = figure; subplot(1,2,1);
    PlotMCC(gca,sweepResults);
    title('MCC')
    AddColorBar()
    subplot(1,2,2);
    PlotMCC(gca,sweepResults,'adjusted');
    title('MCC_a_d_j')
    fH.Position(3:4) = [1000 400];
    AddColorBar()
end

function AddColorBar()
    cb = colorbar('Location','north');
    cb.AxisLocationMode = 'manual';  %cb.YAxisLocation = 'bottom'
    cb.Ticks = cb.Limits(1):0.1:cb.Limits(2);
    cb.Position(1) = cb.Position(1) + cb.Position(3)/3; %x position
    cb.Position(3) = cb.Position(3)/3;  %width
    cb.Position(2) = 1-cb.Position(2) + cb.Position(3);
end


function [bestStart,bestStop] = PlotBestTrainWin(mcc,trainOpt,swOpt,pltColor)
%Find and plot the best train window from Onset, Offset, and Sustained Periods

    % Find best train period
    [mv,mi] = max(mcc(:));
    [ri,ci] = ind2sub(size(mcc),mi);
    bestStart = trainOpt.trainStarts(ri);
    bestStop = trainOpt.trainStops(ci); %relative to Go event


    % Build "best" decoder and train with all data
    m = SetupDecoder(trainOpt.trainFeat,trainOpt.trainInds,trainOpt.trialBounds, trainOpt.startStops, trainOpt.trainLabels, bestStart, bestStop, swOpt.decodeMultistate);
    m.params.enforceTransitionMatrix = 0;
    m.params.smoothDataAmount = [5 0];
    m.params.smoothLikAmount = [5 0];
    m.xVal(); % Get performance numbers
    m.Train(); % Train with all data

    [predState, liklihoods] = m.Test();

    [liklihoods,rawLiks] = m.Predict();

    sampleRate = 50;
 
    startX = find(ismember(trainOpt.trainInds,trainOpt.startStops(:,1))) - 1;
    stopX = find(ismember(trainOpt.trainInds,trainOpt.startStops(:,2))) ;
    line([startX' startX']/sampleRate,[0 1.5],'LineWidth',1.5,'Color',[0.1 0.8 0],'LineStyle',':');
    hold on
    line([stopX' stopX']/sampleRate,[0 1.5],'LineWidth',1.5,'Color',[0.8 0 0],'LineStyle',':');
    
    t = [1:length(liklihoods)]/sampleRate;
    if size(liklihoods,2) > 2 %if showing multistate decoding
        plot(t,liklihoods(:,2:end),'LineWidth',1.5 ); 
    else %binary
        plot(t,liklihoods(:,2),'LineWidth',1.5,'Color','k'); 
    end

    trainWinStart = startX+bestStart;
    trainWinStop = startX+bestStop;
    
    multiColors = lines(8);
    trainLbls = double(trainOpt.trainLabels);

    for i = 1:length(trainWinStart)
        X = [trainWinStart(i) trainWinStart(i) trainWinStop(i) trainWinStop(i)] / sampleRate;
        Y = [1.2 1.3 1.3 1.2];
        if size(liklihoods,2) > 2 %if showing multistate decoding
            pc = multiColors(trainLbls(trainWinStart(i))-1,:);
        else %binary
            pc = pltColor;
        end
        patch(X,Y,pc,'LineStyle','none');
    end
    
    TB = trainOpt.trialBounds;
    TB(1) = []; %remove first one cuz of overlap with yaxis
    lineObjs = line([TB TB]/sampleRate,[0 1.5],'LineWidth',1.5,'Color','k');
    
    ax = gca;
    ax.LineWidth = 1.5;
    ax.FontSize = 13;
    ax.FontWeight = 'bold';

    % put 0's and 4's as xticklabels under appropriate starts and stops
    xticks(1:ax.XLim(2))
    xticklabels('');
    ax.XTickLabels(round(startX'/sampleRate)) = {'Att.'};
    ax.XTickLabels(round(stopX'/sampleRate)) = {'Rel.'};
    ax.XAxis.FontSize = 11;
    ax.XTickLabelRotation = 0;
    
    ax.XAxis.TickLength = ax.XAxis.TickLength/2;
    ax.XAxis.TickDirection = 'none'; %remove xticks
%     ax.XMinorTick = 'on';
%     ax.XAxis.MinorTickValues = 1:ax.XLim(2);
%     ax.XAxis.TickLength(2) = 0;

    %remove 1.5 tick in ylabel
    ax.YTick = [0 1];
    ax.YAxis.TickDirection = 'out';
    ax.YAxis.TickLength = ax.YAxis.TickLength /2;
    
    nTrialsToShow = 5;
    xlim([0 stopX(nTrialsToShow)/sampleRate+2]) %just show 5 example trials

    ylim([0 1.5]);
%     xlabel('sec.')
    ylabel('LDA Likelihood')

    % fix Color Order of likelihood plots in case of subset of gestures used:
    uLbls = unique(trainLbls);
    uLbls(uLbls == 1) = []; % remove no_action label
    liksColorOrder = multiColors(uLbls-1,:);
    ax.ColorOrder = liksColorOrder;
    
end