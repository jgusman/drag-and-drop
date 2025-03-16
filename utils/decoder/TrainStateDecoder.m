

function m = TrainStateDecoder(features, taskInfo, decoderParams, xVal)
% No action is in the labels, select features, then calibrate

    if nargin < 4
        xVal = 1;
    end

    trialBounds = cumsum( [1; diff( taskInfo.startStops,  [], 2 )+1 ] );

    dfParams.kFoldIter = 10;
    dfParams.kFold = 5;
    dfParams.trialBoundaries = trialBounds;
    dfParams.xValVerbosity = 2;
    
    params = MergeBintoA(dfParams, decoderParams);

    params.featInds = taskInfo.rankFeat.selectedInds;

    
    trainInds  = RowColon(taskInfo.startStops);
    flatLabels = GetLabelsForStartStops(taskInfo.labels,taskInfo.startStops);
    % Create decoder instance
    m = dMultistate(params, features(trainInds,:), flatLabels);
        
    % Count num per label
    uLbl = categories(taskInfo.labels);
    for ii = 1:length(uLbl)
        cnt(ii) = sum(ismember(taskInfo.labels,uLbl(ii)));
    end
    
    if xVal
        fprintf('\n\nPerforming %d fold cross-validation on %s decoder...\n', m.params.kFold, taskInfo.labelType)
        mcc = nan;
        if min(cnt) > 6
            try
                % Can possibly fail if xVal doesn't have all classes
                m.xVal();
                mcc = m.performance.mcc;
            end
        end
        
        if isnan(mcc)
            try
                m.params.kFold = 3;
                m.xVal();
            catch
                % If fails, set performance value for Assess_Multistate
                m.performance.mcc = nan;
            end
        end
        fprintf('\n\nAverage performance: %.3f ±%.3f mcc\n\n', m.performance.mcc, m.performance.std.mcc)
    end

    m.Train();

    %% Store some of the training info
    m.info.task = taskInfo;

end
