function PlotLatchProbs(sProbs,taskInfo, predState, latchState, selI, selUpDown, s,TRIALNUM,skipNoAction)
    % Plot gesture and attempt probabilities along with decode outputs of Gesture and Attempt decoders
    
    if ~exist('predState','var')
        predState = [];
    end
    if ~exist('latchState','var')
        latchState = [];
    end
    if ~exist('selI','var') || isempty(selI)
        selI = true(size(taskInfo.startStops,1),1);
    end

    if ~exist('selUpDown','var') || isempty(selUpDown)
        selUpDown = [];
    end
    if ~exist('s','var') || isempty(s)
        s = [];
    end

    extraOnsetOffset = [-15 50];
    
    plotTrialInds = find(selI);
    figWH = [3.5 2.4];
    fh = ddFigHelper.CreateFigure(89898,figWH);
    clf

    axArgs = ddFigHelper.GetDefaultAxArgs();
    axArgs.gap = 0.2;
    axArgs.marg_w = 0.1;
    topAx = 0.7;
    probAxs = (1-topAx)/2;
%     axArgs.proportional_h = [topAx probAxs probAxs ];
    axs = ttight_subplot(length(sProbs)+1,1,axArgs); %3,1,'mergeRC',{1:2, 1});

    lr.ind = 1;
    lr.minMax = [1 length(plotTrialInds)];
    lr.plotFun = @LocalPlotV1;
    lr.wrapBounds = true;

    try
        ind = fh.UserData.lr.ind;
        if ind >= lr.minMax(1) && ind <= lr.minMax
            lr.ind = ind;
        end
    end
    % Just for latch
    startTrl = find(plotTrialInds==TRIALNUM); %choosing trial to plot here!
    if ~isempty(startTrl)
        lr.ind = startTrl;
    end

    data.axs = axs; 
    data.latchState = latchState;
    data.predState = predState;
    data.taskInfo = taskInfo;
    
    data.plotTrialInds = plotTrialInds;
    data.sProbs = sProbs;
    data.extraOnsetOffset = extraOnsetOffset;
    data.attempt = [];
    data.noAttempt = [];
    data.infoStr = InfoString('');
    data.s = s;
    data.skipNoAction = skipNoAction;
    
    fh.UserData = [];
    fh.UserData.lr = lr;
    fh.UserData.data = data;
    
    if ~isempty(selUpDown)
        ud.lastInd = 1;
        ud.ind = 1;
        ud.inds = find(selUpDown);
        ud.minMax = [1 length(ud.inds)];
        fh.UserData.ud = ud;
    end
    
    LeftRightPlot(fh); %allows you to use left and right arrow keys to look through different trials (plotTrialInds)
    lr.plotFun(fh); 
end



%% 
function LocalPlotV1(fh)
%%
        
        figWH = [5.4 2.4]; % [3.5 2.4];
        fh = ddFigHelper.CreateFigure(89898,figWH);
        clf
            
        sProbs = fh.UserData.data.sProbs;
        
        axArgs = ddFigHelper.GetDefaultAxArgs();
        axArgs.gap = 0.1;
        axArgs.marg_w = [0.2 0.04];
        axArgs.marg_h = [0.15 0.05];
        topAx = 0.7;
        probAxs = (1-topAx)/2;
    %     axArgs.proportional_h = [topAx probAxs probAxs ];
        axs = ttight_subplot(length(sProbs)+1,1,axArgs); 
        fh.UserData.data.axs = axs;
        
        data = fh.UserData.data;
        taskInfo = fh.UserData.data.taskInfo;
        
        summaryAx = data.axs(1);
        probAxs = data.axs(2:3);
%
        ti = fh.UserData.lr.ind;
        ti = fh.UserData.data.plotTrialInds(ti);
        
        blk = taskInfo.blockNumber(ti);
%         data.infoStr.String = sprintf('Block %d Epoch %d', blk, ti);
        s = fh.UserData.data.s;
        if ~isempty(s)
            fh.Name = sprintf('[%s] Block %d. Trial %d.', s.date, blk, ti);
        end
        
        % gesture order
        gestIndex = find(ismember({sProbs.name},'Gesture'));
        if gestIndex == 2
            sProbs = flip(sProbs);
            gestIndex = 1;
        end

        %%
        probPos = cat(1,probAxs.Position);
        lowY = min(probPos(:,2));
        highY = max(sum(probPos(:,[2 4]),2));
        
        txtY = mean([highY lowY]);
        txtPos = [0.12, txtY];
        ylbl = YLabelString('Probability', txtPos);
        ylbl.FontSize = 10;
        ylbl.Position(2) = (txtY/2) + 0.05;
        %%
        noDecodeColor = ones(1,3).*0.7;
        skipNoAction = fh.UserData.data.skipNoAction;
        lw = 1;
        lwCued = 2;
        catLbls = categories(taskInfo.labels);
        label = taskInfo.labels(ti);
%         TitleString(char(label))
        lblNum = double(label);
        trlInds = RowColon(taskInfo.startStops(ti,:));
        trl_sec = (trlInds-taskInfo.startStops(ti,1))./50;
        pltInds = RowColon(taskInfo.startStops(ti,:) + data.extraOnsetOffset);
        
        plot_sec = (pltInds - taskInfo.startStops(ti,1))./50;
        trialEnd_sec = taskInfo.trialDurationRounded(ti)./50;
        
        
        unitStr = 'normalized';
        txtX = 0;
        unitStr = 'data';
        txtX = plot_sec(1) - abs(plot_sec(1)*0.2);
            
        yGesture = 1.15;
        
        lwDecodeState = 5;
        
        for axI = 1:length(sProbs)
            ax = probAxs(axI);
            axes(ax);
            cla(ax)
            ax.YLim = [0 1];
            ax.XLim = [plot_sec(1) plot_sec(end)];
            
            [isDecode, plotClrs, probs] = GetDecodedState(sProbs(axI), pltInds);
            threshold = sProbs(axI).threshold;
            nDIm = size(probs,2);
            clrs = sProbs(axI).colors;
            
            for ii = 1:nDIm
                if skipNoAction && ii == 1
                    continue
                end
                c = clrs(ii,:);
%                 if ii == lblNum && length(catLbls) == nDIm
%                     pltLw = lwCued;
%                 else
                    pltLw = lw;
%                 end
                
                plot(plot_sec,probs(:,ii),'-','LineWidth',pltLw,'Color',c)
                isOverThresh = probs(:,ii) > threshold;
                plotVal = nan(size(isOverThresh));
                plotVal(isOverThresh) = probs(isOverThresh,ii);
                plot(plot_sec,plotVal,'-','LineWidth',2,'Color',c)
            end

            clrI = find(ismember(ddHelper.decoder.names,sProbs(axI).name));
            c = ddHelper.decoder.colors(clrI,:);
            ax.YLabel.String = sprintf('%s', sProbs(axI).name);
            ax.YLabel.Color = c;
            

            yl = yline(threshold,'r:');

            yl.LineWidth = 0.5;
            xline(0,'k');
            xline(trialEnd_sec,'--k');
        end
%
        txtFs = 9;
        txtY = -0.4;
        
        th = text(0,txtY,'Attempt onset', 'HorizontalAlignment','center', 'Clipping', 'off', ...
            'FontSize', txtFs);
        
        th = text(trialEnd_sec,txtY,'End', 'HorizontalAlignment','center', 'Clipping', 'off', ...
            'FontSize', txtFs);
        
        xlabel('Time (s)')
            
            %%
            
            
            ax = summaryAx;
            ax.XRuler.TickLength(1) = 0.007;
            axes(ax);
            cla(ax)
            ax.YColor = 'none';
            ax.YLimMode = 'auto';
            ax.XColor = 'none';
            ax.XLim = [plot_sec(1) plot_sec(end)];
            xline(0,'k');
            xline(trialEnd_sec,'--k');
            
            labelH = 0.12;
            trialLblLW = 15;
            stateLW = 4;
            
            yTop = 1;

            
            ys = linspace(0.1,0.8,3);
            yAttemptDecode = ys(1); 
            yGestureDecode = ys(2);
            yOut = ys(3);
            yBottom = 0;
            
            
            decodeH = 0.1;
            noDecodeColor = ones(1,3).*0.7;
            for ii = 1:length(sProbs)
                [isDecode, plotClrs] = GetDecodedState(sProbs(ii), pltInds);
                
                if strcmpi(sProbs(ii).name,'Gesture')
                    yVal = yGestureDecode;
                else
                    yVal = yAttemptDecode;
                end
                plotVal = repmat(yVal,length(isDecode),1);

                plotClrs(~isDecode,:) = repmat(noDecodeColor, sum(~isDecode),1);
                h = PlotColorLine(plot_sec, plotVal,plotClrs, 'LineWidth',stateLW, 'Clipping','off');
            end

            % ----------------
            % Latch
            % ----------------
            isLatched = data.latchState(pltInds);
            latchOnsetOffset = diff([false; isLatched]);
            latchOnsets = find(latchOnsetOffset>0);
            latchClr = [0.1919, 0.6427, 0.3419];
            for ii = 1:length(latchOnsets)
                latchX = plot_sec(latchOnsets(ii));
                xl = xline(latchX,'-'); %,'Latch enabled')
                xl.Color = latchClr;

                
                % Centered text

                plot(latchX,yTop,'v', 'MarkerFaceColor', latchClr, 'MarkerEdgeColor','none');
                latchTxtX = latchX - 0.05;
                latchTxtY = yTop + 0.18;
                text(latchTxtX,latchTxtY, 'Latch enabled', 'HorizontalAlignment', 'left', 'FontSize', 7, 'Clipping','off')

                
            end
            
            latchOffsets = find(latchOnsetOffset<0);
            offsetClr = [1 0 0];
            for ii = 1:length(latchOffsets)
                latchX = plot_sec(latchOffsets(ii)-1);
                xl = xline(latchX,'-'); %,'Latch disabled')
                xl.Color = offsetClr;
                plot(latchX,yTop,'v', 'MarkerFaceColor', offsetClr, 'MarkerEdgeColor','none');
                latchTxtX = latchX - 0.05;
                latchTxtY = yTop + 0.18;
                text(latchTxtX,latchTxtY, 'Latch disabled', 'HorizontalAlignment', 'left', 'FontSize', 7, 'Clipping','off')
            end
            

            % ----------------
            % Output
            % ----------------
            outputState = data.predState(pltInds);
            isDecode = outputState>1;
            plotVal = repmat(yOut,size(outputState));
            plotClrs = sProbs(gestIndex).colors(outputState,:);

            % Set no decode to be grey
            plotClrs(~isDecode,:) = repmat(noDecodeColor, sum(~isDecode),1);
            h = PlotColorLine(plot_sec, plotVal,plotClrs, 'LineWidth',stateLW, 'Clipping','off');
            

            % ----------------
            % Text
            % ----------------
            txtX = -0.5;
            txtStateY = 1.1;
            th = text(txtX-0.3,txtStateY,sprintf('Decoded state'),'Units',unitStr,'HorizontalAlignment','center','FontWeight', 'bold','Clipping','off');
            
            
            txtStrs = {'Attempt', 'Gesture', 'Latch'};
            for ii = 1:length(txtStrs)
                clrI = find(ismember(ddHelper.decoder.names,txtStrs{ii}));
                c = ddHelper.decoder.colors(clrI,:);
                txtY = ys(ii);
                th = text(txtX,txtY,txtStrs{ii},'Units',unitStr,'HorizontalAlignment','right','Color', c);
            end
%             th = text(txtX,yGestureDecode,'Gesture','Units',unitStr,'HorizontalAlignment','right');
%             th = text(txtX,yAttemptDecode,'Attempt','Units',unitStr,'HorizontalAlignment','right');
%             th = text(txtX,yOut,'Decoder','Units',unitStr,'HorizontalAlignment','right');

%             th = text(txtX,yLatch,'Latch enable','Units',unitStr,'HorizontalAlignment','right');
            
%             th = text(txtX,yGestError,'Gesture error','Units',unitStr,'HorizontalAlignment','right');
%             th = text(txtX,yOutError,'Output error','Units',unitStr,'HorizontalAlignment','right');
            
            plot(0,yTop);
            plot(0,yBottom);
        
end

function [isValid, plotClrs, probs,decodeState] = GetDecodedState(sProb, pltInds)
        probs = sProb.probs(pltInds,:);
        threshold = sProb.threshold;

        [mv,decodeState] = max(probs,[],2);
        decodeState(mv<threshold) = 1;
        isValid = decodeState>1;

        plotClrs = sProb.colors(decodeState,:);
end
