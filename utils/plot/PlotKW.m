
function ax = PlotKW(kwPerDurationPerSesAllFeat,sigAlpha,figNum,figWH, clrs,separateFeatsTypes,uDurations)
% Plots average significant features (averaged across sessions) for each hold duration.


if nargin < 6
    separateFeatsTypes = 0; % if we want to plot feat types separately
end

if nargin < 7
    uDurations = [50 100 200];
end

uDurations_sec = uDurations./50;
kws = kwPerDurationPerSesAllFeat;

nSessions = length(kws);
nDurations = length(kws{1});


nFeat = size(kws{1}{1}.p,2);

if separateFeatsTypes
    numFeatTypes = round(nFeat/192);
    featPerType = 192;
    lineWidth = 1.5;
else
    numFeatTypes = 1;
    featPerType = nFeat;
    lineWidth = 2;
end


% Set up figure
ddFigHelper.CreateFigure(figNum,figWH)
clf
axArgs = ddFigHelper.GetDefaultAxArgs();
axArgs.marg_w = [0.13 0.05]; % Update axes width boundaries
axArgs.marg_h = [0.2 0.08]; % Update axes height boundaries
axArgs.xTickOffset = -2; % Make x ticks closer to axes, I think.

ax = ttight_subplot(1,1,axArgs);


hold on

% Plot the average percent selective for each duration
% Averaged across all sessions tested
for di = nDurations:-1:1   %plot them so that 1 sec and 2 sec are in front fo 4 sec for clarity
    val = [];

    for fT = 1:numFeatTypes

    featInds = ((fT-1)*featPerType+1) : (fT*featPerType);

    for jj = 1:nSessions
        kwPlot = kws{jj}{di};
        % Percent significant across all features tested
        % You may want to break this into subsets of features
        val(:,jj) = mean(kwPlot.p(:,featInds)<sigAlpha,2)*100;
    end
    % Average across sessions
    val = mean(val,2); 
    
    % some mild smoothing to immprove plot readability
    val = movmean(val,5);     

    slideTime = kwPlot.slideWinInds./50;
    
    plotI = (di-1)*numFeatTypes + fT;

    plot(slideTime,val,'LineWidth',lineWidth, 'Color', clrs(plotI,:))

    dur = uDurations(di);
    [peakVal,peakInd] = max(val);
    fprintf('(%d Sec Trials) Peak: %.2f%% SelFeats at %.2f sec. End: %.2f%% at %.2f sec.\n ',dur/50,peakVal,slideTime(peakInd),val(slideTime==dur/50),dur/50)


    valMaxs(di) = max(val);
    end
end


% Now plot onset and offsets lines
drawnow;
if max(valMaxs) + 1 > ax.YLim(2)
    ax.YLim(2) = ax.YLim(2) + 5;
end
% ax.YLim(2) = 30;
orgYLim = ax.YLim;

% Plot black line at time 0 (attempt onset)
% Not using xline because it is always in front of plotted lines
% xl = xline(0,'k');

xs = repmat(0,1,2);
ys = orgYLim; %[0 ax.YLim(2)];
xl = plot(xs,ys,'k','LineWidth', 0.5, 'Clipping', 'off');
uistack(xl,'bottom')


% Plot attempt offset line with triangle for each duration.
for di = 1:length(uDurations)
    xs = repmat(uDurations_sec(di),1,2);

    xl = plot(xs,ys,'-','color',clrs(di,:),'LineWidth', 0.15, 'Clipping', 'off');
    uistack(xl,'bottom')
    xl = plot(xs(end),ys(end),'v','color',clrs(di,:), 'MarkerFaceColor',clrs(di,:),'LineWidth', 1, 'MarkerEdgeColor','none', 'Clipping', 'off');
    uistack(xl,'bottom')
    
%     xl = xline(uDurations_sec(di),'-');
%     xl.Color = clrs(di,:);
end


% Force y limits and x ticks to be on second boundaries 
ax.YLim = orgYLim;
ax.XTick = RowColon(ax.XLim);

xlabel('Trial Time (sec)','FontName','Arial','FontSize',10) 
ylabel('Selective Features (%)','FontName','Arial','FontSize',7.5)
end

