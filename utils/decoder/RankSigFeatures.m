function [rankFeat, selectedInds] = RankSigFeatures(features, labels, trStStp, excludeFeatInds, sigAlpha, minMax, selTrialsKW)
% rankFeat = RankSigFeatures(features, taskInfo, excludeFeatInds, sigAlpha)
% Rank features by using the Minimum Redundancy Maximum Relevance algorithm

if nargin < 4
    excludeFeatInds = [];
end
if nargin < 5
    sigAlpha = 0.01;
end
if nargin < 6
    minMax = [];
end
if nargin < 7
    selTrialsKW = true(size(trStStp,1),1); 
end
	
% Subselect significant features before computing 
uLabels = unique(labels);
if length(uLabels) <=2 && any(ismember(uLabels,'no_action'))
    ignoreNoActionForFeatureRanking = false;
else
    ignoreNoActionForFeatureRanking = true;
end

% Subselect significant features before computing 
nSteps = size(features,1);

% Remove no_action trials when finding selective features
if ignoreNoActionForFeatureRanking
    if size(labels,1) == nSteps
        toRmv = any(ismember( labels(trStStp),'no_action'),2);
        trStStp(toRmv,:) = [];
        selTrialsKW(toRmv) = [];
    else
        toRmv = any(ismember( labels,'no_action'),2);
        labels(toRmv,:) = [];
        trStStp(toRmv,:) = [];
        selTrialsKW(toRmv) = [];
    end
end

trainInds = RowColon(trStStp);

if size(labels,1) == nSteps
    labelPerTrial = labels(trStStp(:,1));
    labelPerStep  = labels(trainInds);
else
    labelPerTrial = labels;
    labelPerStep  = GetLabelsForStartStops(labels, trStStp);
end

    labelsKW = labelPerTrial(selTrialsKW); % e.g. just want the isAction trials for gesture decoder of latch
    trStStpKW = trStStp(selTrialsKW,:);
    kwp     = KWperFeature(features, labelsKW, trStStpKW, excludeFeatInds);
sigInds = find(kwp<sigAlpha);



% Run Minimum Redundancy Maximum Relevance algorithm

if ~isempty(sigInds)
    nSig = length(sigInds);
    if ~isempty(minMax) && (nSig > max(minMax))
        [idx,scores] = fscmrmr(features(trainInds,sigInds),labelPerStep);
    else
        [scores, idx] = sort(kwp(sigInds));
    end
else
    [scores, idx] = sort(kwp(sigInds));
end


%%

%%
% Combine for output: ranked feat inds, score, and kruskal wallis p values 
scoreVals = nan(size(features,2),1);
scoreVals(sigInds) = scores;

rankFeat.inds  = sigInds(idx);
rankFeat.score = scoreVals;
rankFeat.kwp   = kwp;

selectedInds = GetSelectedInds(rankFeat, minMax);
rankFeat.selectedInds = selectedInds;


%% Debug plot
debugPlot = 0;
if debugPlot
    DebugPlot(features,rankFeat, labelPerTrial, trStStp);
end
end

function selectedInds = GetSelectedInds(rankFeat, minMaxFeat)
selectedInds = rankFeat.inds;
nSigFeat = length(selectedInds);


if ~isempty(minMaxFeat)
    if nSigFeat < minMaxFeat(1)
        nFeatNeeded = minMaxFeat(1)-nSigFeat;
        [~,si] = sort(rankFeat.kwp);
        unusedSortedFeatInds = setdiff(si,selectedInds, 'stable');
        newFeat = unusedSortedFeatInds(1:nFeatNeeded);
        selectedInds = cat(1, selectedInds(:), newFeat);
    elseif nSigFeat > minMaxFeat(2)
        selectedInds = selectedInds(1:minMaxFeat(2));
    end
end
selectedInds = sort(selectedInds);
end

function DebugPlot(features,rankFeat, labelPerTrial, startStops)
%     figure
%     bar(rankFeat.inds, rankFeat.score(rankFeat.inds))
%     ax = gca;
%     ax.YScale = 'log';
%     xlabel('Predictor rank')
%     ylabel('Predictor importance score')
%     title(sprintf('%s feature predictor importance score', taskInfo.labelType))
%     
    %%
    plotWindow = -25:75;
    plotFeatRank = [1 2 10 20 30 40 nFeat];
    nSigFeat = length(rankFeat.inds);
    isSmallSigFeat = plotFeatRank> nSigFeat;
    if any(isSmallSigFeat)
        plotFeatRank(isSmallSigFeat) = [];
        plotFeatRank(end+1) = nSigFeat;
    end
    
    plotFeat = rankFeat.inds(plotFeatRank);
    smthFeat = GaussSmooth2(features(:,plotFeat), 80);
    nSp = length(plotFeatRank);
    
    figure(100)
    clf
    tsp = tight_subplot(3,2);
    nSp = min(length(tsp),nSp);
    for ii = 1:nSp
        axes(tsp(ii));
%         subplot(nSp,1,ii)

        PlotPeristimFeature(smthFeat(:,ii),labelPerTrial, startStops, plotWindow)
        title(sprintf('Feature %d (Rank %d)', plotFeat(ii), plotFeatRank(ii)))
    end
    InfoString(sprintf('Ranked by: %s', type));
    drawnow

end

function PlotPeristimFeature(feature,labelPerTrial, startStops, plotWindow)
uLbls = unique(labelPerTrial);
for ii = 1:length(uLbls)
    selI = ismember(labelPerTrial,uLbls(ii));
    fEpoch = GetEpochOfData(feature, [],startStops(selI,1), plotWindow,1);
    [~,hh(ii)] = shadedErrorBar(plotWindow./50, fEpoch', 'norm', 'linewidth', 2);
end
legend(hh,uLbls)

end

function kwp = KWperFeature(features, labelPerTrial, trStStp, excludeFeatures)

%% Get trials to apply Kruskal-Wallis test

nSteps = size(features,1);


%% Average features for each trial start stop
trlLen = diff(trStStp,[],2);
if any(trlLen<1)
    error('Found %d trials that were < 1 step in length', sum(trlLen<1))
end

nTrials = size(trStStp,1);
nFeat = size(features,2);
avgData = nan(nTrials,nFeat);

for ii = 1:size( trStStp, 1)
    inds = trStStp(ii,1):trStStp(ii,2);
    avgData(ii,:) = mean(features(inds,:));
end


%% Compute Kruskal-Wallis p value per feature
skipFeat = ismember(1:nFeat,excludeFeatures);

kwp = ones(nFeat,1);
for f = 1:nFeat
    if skipFeat(f)
        continue
    end
    d = avgData(:,f);
    kwp(f) = kruskalwallis(d(:), labelPerTrial(:), 'off');
end


end


