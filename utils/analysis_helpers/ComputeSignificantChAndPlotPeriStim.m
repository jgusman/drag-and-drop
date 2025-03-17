 function [kwP,kwOpt,handles] = ComputeSignificantChAndPlotPeriStim(features, varargin)
% [kwP,kwOpt,handles] = ComputeSignificantChAndPlotPeriStim(features, labels, events, winInds, slideWinInds, baselineEvents, baselineStatWinInds)
% where the following are optional: winInds, slideWinInds, baselineEvents, baselineStatWinInds
% 
% Or pass in as a struct
% 
% [kwP,kwOpt] = ComputeSignificantChAndPlotPeriStim(features, kwOpt)
%
% winInds       relative inds around selectInds used to compute KW per slide step
%               For example, 0:14 will be the first 15 steps after each
%               selectInd.

defaultOpt.events = [];
defaultOpt.labels = [];
defaultOpt.winInds = -14:0;
defaultOpt.slideWinInds = [];
defaultOpt.baselineEvents = [];
defaultOpt.baselineStatWinInds = -14:0;
defaultOpt.skipStepsPerCh = [];
defaultOpt.plotOn = true;
defaultOpt.plotWinInds = -50:50;
defaultOpt.plotTestAlpha = 0.01;
defaultOpt.plotColors = []; % Should be [nConditions x 3]

%% Get default params if called with no inputs.
if nargin == 0
    [kwP,kwOpt] = deal(defaultOpt);
    handles = [];
    return;
end


%% Handle inputs


nArgs = length(varargin);
if isstruct(varargin{1})
    kwOpt = varargin{1};
else
    argumentsInOrder = {'labels', 'events',  'winInds', 'slideWinInds', 'baselineEvents', 'baselineStatWinInds', 'skipStepsPerCh'};
    for ii = 1:nArgs
        kwOpt.(argumentsInOrder{ii}) = varargin{ii};
    end
end

kwOpt = MergeBintoA(defaultOpt,kwOpt);


if isempty(kwOpt.slideWinInds)
    kwOpt.slideWinInds = 0;
end

if isfield(kwOpt, 'statWinInds')
    error('Using please use ''winInds'' instead of ''statWinInds''.');
end

%% Run analysis


if isvector(kwOpt.baselineEvents) || isempty(kwOpt.baselineEvents)
    kwP = ComputeSigCh(features, kwOpt.labels, kwOpt.events, kwOpt.winInds, kwOpt.slideWinInds, kwOpt.baselineEvents, kwOpt.baselineStatWinInds, kwOpt.skipStepsPerCh);
else
    kwP = kwOpt.baselineEvents;
    kwOpt.baselineEvents = [];
end



%%

if kwOpt.plotOn
    sFeat = GaussSmooth2(features,80);
    minMaxSlideInds = [min(kwOpt.slideWinInds) max(kwOpt.slideWinInds)];
    minMaxSlideInds = [minMaxSlideInds(1) + min(kwOpt.winInds) minMaxSlideInds(2) + max(kwOpt.winInds)];
    
    plotWinInds = min(minMaxSlideInds(1),min(kwOpt.plotWinInds)):max(minMaxSlideInds(2), max(kwOpt.plotWinInds));
%     plotWinInds = -25:75;
    

    
    
    if size(kwP,1) > 1
        % Plot multiple time steps
%         sigChPerStep = sum(kwP<kwOpt.plotTestAlpha,2);
%         figure
%         plot(slideWinInds+max(statWinInds), sigChPerStep)
%         xlabel('Steps')
%         ylabel('Num sig features')
%         title(sprintf('Sig features per step (p<%.4f)', kwOpt.plotTestAlpha))
    end
    
    
    
    handles = PlotSigFeat(sFeat, removecats(kwOpt.labels), kwOpt.events, plotWinInds, kwP, kwOpt);
    drawnow
else
    handles = [];
end

end

function kwP = ComputeSigCh(features, labels, selectInds, statWinInds, slideWinInds, baselineEvents, baselineStatWinInds, skipStepsPerCh)


avgPerTrial = 1;

isLprintfAvailable = ~isempty(which('lprintf'));
if isLprintfAvailable
    printFun = @lprintf;
else
    printFun = @fprintf;
end

%% Skip steps
if ~isempty(skipStepsPerCh)
    features(skipStepsPerCh) = nan;
    
end

%% Baseline

b = [];
baselineLabel = categorical();
isBaseline = ~isempty(baselineEvents);
if isBaseline
    baseND = GetEpochOfData(features, [], baselineEvents, baselineStatWinInds, 1);
    baselineLabel = categorical( repmat({'Baseline'}, length(baselineEvents),1) );
    
    if avgPerTrial
        baseND = squeeze(nanmean(baseND,1));
    end
else
    baseND = [];
end

%%
if iscategorical(labels)
    labels = removecats(labels);
end

%%
nSweep = length(slideWinInds);
nDim = size(features,2);
kwP = ones(nSweep, nDim);

if isLprintfAvailable
    printFun();
end
for swpI = 1:nSweep
    if nSweep > 1
        printFun('calculated %d out of %d windows\n', swpI, nSweep);
    end
    sweepEventInds = selectInds + slideWinInds(swpI);
    [statTestND, lbl] = GetEpochOfData(features, labels, sweepEventInds,statWinInds ,1);
    if avgPerTrial
        
        statTestND = nanmean(statTestND,1);
        lbl = lbl(1,:);
    end

    for ii = 1:nDim
        
        d = statTestND(:,:,ii);
        if isBaseline
            b = baseND(:,ii);
        else
            b = [];
        end
        

        kwP(swpI,ii) = kruskalwallis( cat(1,d(:),b(:)), ...
                                      cat(1,lbl(:),baselineLabel(:)), 'off' );
    end

end
end







function handles = PlotSigFeat(smthFeat, cuePerStep, selectInds, slideWinInds, kwP, kwOpt)
%% Plot
usedBaseline = ~isempty(kwOpt.baselineEvents);

isSig = kwP<kwOpt.plotTestAlpha;
isMultipointTest = size(kwP,1) > 1;

if isMultipointTest
    sigChPerStep = sum(isSig);
    [sigChPerStep,si] = sort(sigChPerStep, 'descend');
    sigInds = si(sigChPerStep>0);
else
    sigInds = find(isSig);
    [~,si] = sort(kwP(isSig));
    sigInds = sigInds(si);
end



if usedBaseline
    [baselineND] = GetEpochOfData(smthFeat, [], kwOpt.baselineEvents, kwOpt.baselineStatWinInds,1);
end

[epochND, epochLabel, epoch] = GetEpochOfData(smthFeat, cuePerStep, selectInds, slideWinInds,1);

%%

firstTrInd = find(slideWinInds==0); if isempty(firstTrInd), firstTrInd = 1; end
labelPerTrial = cellstr( epochLabel(firstTrInd, : ) );



maxPlot = 9;
nPlot = length(sigInds);
if nPlot == 0
    handles = [];
    return;
end

nSubPlot = min( ceil(sqrt(nPlot)), sqrt(maxPlot) );
% Set nPlot to be the max number of values we plot
if nSubPlot^2 < nPlot
    nPlot = nSubPlot^2;
end
% plot in order
sortSigInds = sort(sigInds(1:nPlot));
    
pltTime = (slideWinInds)./50;
testTime = kwOpt.slideWinInds./50; % What points were tested
if isscalar(testTime)
    testTime = testTime + kwOpt.winInds./50;
end

[~, uLbl] = grp2idx(cuePerStep);
contClrs = cmap(length(uLbl));
if isequal(size(kwOpt.plotColors),size(contClrs)) && size(kwOpt.plotColors,2)==3
    contClrs = kwOpt.plotColors;
end
%%    
fh = MaximizeFig(109,2);
clf
for ii = 1:nPlot
    ax = subplot(nSubPlot,nSubPlot,ii);
    sp(ii) = ax;
    ch = sortSigInds(ii);
    clear hh
    for ui = 1:length( uLbl )

        selI = ismember(labelPerTrial, uLbl(ui));
        nd = epochND(:,selI,ch);
        
        clr = contClrs(ui,:);
        [~,hh(ui)] = shadedErrorBar(pltTime, nd', 'norm', 'linewidth', 2, 'color', clr );
%         hh(ui) = plot(pltTime, mean(nd), 'LineWidth', 2, 'color', clr );
        
    end
    
    if usedBaseline
        %%
%         baselineClr = ones(1,3).*0.4;
        d = baselineND(:,:,ch);
        [u, ~,ci] = normfit(d(:));
        bh = plot(pltTime([1 end]),[ci(1) ci(1)], '--k');
        plot(pltTime([1 end]),[ci(2) ci(2)], '--k');
%         ci = u - ci(1);
%         [~,bh] = shadedErrorBar(pltTime([1 end]),[u u], [ci ci], 'linewidth', 1, 'color', baselineClr );
        
    end
    
    % If multiple time steps were tested, plot where significant.
    if isMultipointTest
        ph = PlotSigEpoch(ax, testTime, ch, kwP, kwOpt);
    end
    
    stStpColor = ones(1,3).*0.65;
    xl = xline(testTime(1),'-', 'start test');
    xl.Color = stStpColor;
    xl = xline(testTime(end),'-', 'stop test');
    xl.Color = stStpColor;
    
    
    xlim([min(pltTime) max(pltTime)]);
    
    titleStr = sprintf('Feat. %d', ch);
    title(titleStr)
    
    if ii == 1
        legH = hh;
        legStr = uLbl;
        if usedBaseline
            legH = [legH bh];
            legStr{end+1} = 'Baseline';
        end
        if isMultipointTest
            legH = [legH ph];
            legStr{end+1} = 'Significant points';
        end
        lh = legend(legH, legStr);
        lh.Interpreter = 'none';
    end
    if nPlot == 9
        if ii == 4
            ylabel('Feature value')
        end
        if ii == 8
            xlabel('Seconds')
        end
    end
    
    
end
infStr = InfoString(sprintf('%d features p < %.1e', length(sigInds), kwOpt.plotTestAlpha));

lh.Position(1:2) = [0.0040    0.9 - lh.Position(4)];
handles.fig = fh;
handles.ax = sp;
handles.infoStr = infStr;

% legend([hh], uLbl);

end

function ph = PlotSigEpoch(ax, testTime, ch, kwP, kwOpt)
    sigClr = [0.5006, 0.8133, 0.6423];
    
    lw = 3;
    ms = 3;

    isSig = kwP(:,ch) < kwOpt.plotTestAlpha;
    sigY = ax.YLim(2) - ax.YLim(2)*0.05;
    pltSig = nan(size(testTime));
    pltSig(isSig) = sigY;
    

%     fprintf('%d is sig %d steps at Y=%.2f\n', ch, nSigPoints, sigY);
    %%
    % Find isolated points
    [ws,wststp] = FindConsecutiveOnes(isSig,1);
    isolatedPoints = find(ws == 1);
    if ~isempty(isolatedPoints)
        % Plot square markers for single points since these do not show up
        % in a nan line
        pltI = wststp(isolatedPoints,1);
        ph = plot(testTime(pltI), pltSig(pltI), 's', 'Color', sigClr, 'LineWidth',lw, 'MarkerFaceColor', sigClr, 'MarkerSize', ms);
    end
    ph = plot(testTime, pltSig, '-', 'Color', sigClr, 'LineWidth',lw, 'MarkerFaceColor', sigClr, 'MarkerSize', ms);
%%
    
end