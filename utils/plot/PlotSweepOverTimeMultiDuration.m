function ax = PlotSweepOverTimeMultiDuration(swp,axs,plotOpt)
    
    if ~exist('axs','var') || isempty(axs)
        ddFigHelper.CreateFigure(32239)
        clf
        axArgs = ddFigHelper.GetDefaultAxArgs();
        axArgs.marg_w = 0.1;
        axArgs.marg_h = 0.12;
        axs = ttight_subplot(1,1,axArgs);
    end

    dfPlotOpt.showMcc = 0;
    dfPlotOpt.forceYLim = [];
    dfPlotOpt.mainM = '-';
    dfPlotOpt.showCI = 1;
    dfPlotOpt.showLegend = 1;
    dfPlotOpt.lw = 2;
    dfPlotOpt.lwPerSes = 0.5;
    dfPlotOpt.plotPerSession = 0;
    dfPlotOpt.perSessionMarkers = {'--',':'}; % if plotPerSession is true. A linestyle/marker for each session
    dfPlotOpt.clrs = [0.2200    0.6364    0.9892
                      0.1050    0.4694    0.7781
                      0.0600    0.2682    0.4446];

    if ~exist('plotOpt','var')
        plotOpt = [];
    end

    plotOpt = MergeBintoA(dfPlotOpt,plotOpt);
    clrs = plotOpt.clrs;


%     m = {'--',':'};
% 
%     m = m(1:length(swp));


    %%
    % swp = FixOffset(swp);
    allOffsets = swp(1).offsets;
    uTrialDuration = swp(1).durations;
    nDurs = length(uTrialDuration);
    nFeatTypes = size(swp,1);
    stepSize = swp(1).info.stepSize;

    axes(axs(1));

    for di = nDurs:-1:1 %plot them so that 1 sec and 2 sec are in front fo 4 sec for clarity
        dur = uTrialDuration(di);
        
        for fT = 1:nFeatTypes  % plot feature types separately if given many

            allPerfs = cat(1,swp(fT,:).perf);
            %     allPerfs = swp(swpI).perf;
            perfs = cat(1,allPerfs{:,di});
            vals = [];
            for ii = 1:size(perfs,1)
                if ~plotOpt.showMcc
                    for jj = 1:size(perfs,2)
%                         vals(ii,jj) = mean(cat(2,perfs(ii,jj).recall(2:end)),1).*100;
                        vals(ii,jj) = mean(cat(2,perfs(ii,jj).recall),1).*100;
                    end
                else
                    vals(ii,:) = (cat(2,perfs(ii,:).mcc));
                end
            end
    
            offsetT = allOffsets{di}./50;
            avgVals = mean(vals,1);

            plotI = (di-1)*nFeatTypes + fT;
            c = clrs(plotI,:);
            hh(plotI) = plot(offsetT,avgVals,plotOpt.mainM,'LineWidth', plotOpt.lw, 'color', c);
    
            [peakVal,peakInd] = max(avgVals);

            if isempty(avgVals(offsetT==dur/50))
                [val,idx]=min(abs(offsetT-dur/50));
                endTime=offsetT(idx);
            else
                endTime = dur/50;
            end

            fprintf('(%d Sec Trials) Peak: %.2f%% Correct at %.2f sec. End: %.2f%% at %.2f sec.\n ',dur/50,peakVal,offsetT(peakInd),avgVals(offsetT==endTime),endTime)
    
            if plotOpt.plotPerSession
                cSess = min(c.*1.1,1);
                for ii = 1:size(perfs,1)
                    plot(offsetT,vals(ii,:),plotOpt.perSessionMarkers{ii},'LineWidth', plotOpt.lwPerSes, 'color', cSess)
                end
            end

        end

    end


    if ~plotOpt.showMcc && ~isempty(plotOpt.showCI) && ~isequal(plotOpt.showCI,0)
        
        for fT = 1:nFeatTypes
            swpInfo = [swp(fT,:).info];
            lbls_cell = cat(1, swpInfo.labels);
            lbls = cat(1,lbls_cell{:});
            chanceCI(fT,:) = ComputeChanceOnAttemptOrGestures(swp,lbls);
            % chanceVal = 1 / (length(swp(1).decoder.info.train.uStates)-1)*100;
        end
        chanceCI = mean(chanceCI,1);

        if ischar(plotOpt.showCI)
            errMarker = plotOpt.showCI;
        else
            errMarker = '--';
        end
        yline(chanceCI(1),errMarker);
        yline(chanceCI(2),errMarker);
    end

    xlabel('Trial Time (sec)')
    if plotOpt.showMcc
        ylabel('MCC')
    else
%         ylabel(sprintf('Sensitivity (%%)'))
        ylabel('Percent Correct')
    end

    axis tight
    drawnow;
    ax = gca;
    if ~isempty(plotOpt.forceYLim)
        ax.YLim = plotOpt.forceYLim;
        orgYLim = plotOpt.forceYLim;
    else
        orgYLim = ax.YLim;
    end
    
    % Draw Onset Line
    xs = repmat(0,1,2);
    ys = orgYLim; %[0 ax.YLim(2)];
    xl = plot(xs,ys,'k','LineWidth', 0.5, 'Clipping', 'off');
    uistack(xl,'bottom');
    % xline(0,'k');

    % Draw Trial End Lines
    for di = 1:length(uTrialDuration)
        xs = repmat(uTrialDuration(di)./50,1,2);

        xl = plot(xs,ys,'-','color',clrs(di,:),'LineWidth', 0.15, 'Clipping', 'off');
        uistack(xl,'bottom');
        xl = plot(xs(end),ys(end),'v','color',clrs(di,:), 'MarkerFaceColor',clrs(di,:),'LineWidth', 0.5, 'MarkerEdgeColor','none', 'Clipping', 'off');
        uistack(xl,'bottom');
        %     xl = xline(uTrialDuration(di)./50,'-');
        %     xl.Color = clrs(di,:);
    end
    ax.YLim = orgYLim;

    if plotOpt.showLegend
        legStr = sprintfc('%.0f sec hold', uTrialDuration./50);
        if size(swp,2) > 1 && plotOpt.plotPerSession
            for ii = 1:size(swp,2)
                hh(end+1) = plot(nan,[plotOpt.perSessionMarkers{ii} 'k']);
                legStr{end+1} = swp(1,ii).info.s.date;
            end
        end
        hh(end+1) = plot(nan,'-k');
        legStr{end+1} = 'Attempt onset';

        hh(end+1) = plot(nan,'--k');
        legStr{end+1} = 'Chance'; % '95 CI';
        lh = legend(hh,legStr);
        lh.AutoUpdate = 'off';
    end

end


function chanceCI = ComputeChanceOnAttemptOrGestures(swp,lbls)
    ciAlpha = 5;
    N = 10e3;
    ulbls = unique(lbls);
    if length(ulbls) > 2
        % Gestures
        lbls(ismember(lbls,'no_action')) = [];
        chanceCI = GetChanceCI(lbls, ciAlpha, N);
    else
        % 2 class (attempt only)

        % Remove no_action labels since we only care about attempt trials
        lbls(ismember(lbls,'no_action')) = [];
        lbls = lbls(:);

        shuffledLabels = lbls;
        nLblTrials = length(shuffledLabels);
        prctScale = 100/nLblTrials;
        chanceCorrect = nan(N,1);

        for iter = 1:N
            shuffledLabels = lbls;
            % Randomly set labels to be no_action
            setLbls = randi(2,nLblTrials,1) == 1;
            shuffledLabels(setLbls) = 'no_action';
            chanceCorrect(iter) = sum(lbls == shuffledLabels)*prctScale;
        end

        if isscalar(ciAlpha)
            prctileVals = [ciAlpha/2 100-ciAlpha/2];
        else
            prctileVals = ciAlpha;
        end
        chanceCI = prctile(chanceCorrect, prctileVals);
    end
end
