classdef dMultistate < Decoder
    % This is a multistate hmm discrete decoder
    % Run 'help dMultistate' to get a summary of all functions in class and 
    % parameters 
    
    % ToDo:
    % Add feature selection option using rankfeatures
    %
    % By depending on another decoder class for the HMM, we may want to
    % make the coefficients more static and only allow the user to set them
    % with specific set functions to prevent a mismatch of coefficient
    % values.
    
    properties
        data  % neural data, [nSteps x nFeat]
        labels % discrete labels for each class,  [nSteps x 1]
        
        eventInds % eventInds.train, eventInds.test. Default to empty (use all).
        performance % struct containing performance results
        coefficients % struct with parameters, coefficients for the decoder
        params % Defines how to calc the multistate decoder        
        info % Saved predicted data used for plotting in Summarize function 
        stateNames 
    end
    
    
    
    methods
        function obj = dMultistate( params, data, labels )
            % params can be user specified params to overwrite any default params.
            % Or params can be structs from a session block.
            %   slc struct with .sSLC and .task
            %   singleBlock struct
            %   X (not yet supported) sFilter struct 
            % 
            % All inputs are optional (but data and labels must be filled
            % out eventually)
        
            % Default parameters
            defaultParams.dimReductionMethod           = 'lda'; % lda, pca, none
            defaultParams.ldaRegularization            = 0.5;
            defaultParams.featInds                     = []; % What features from data (2nd dim) to use. If empty, all are used.
            defaultParams.selectFeaturesMinMax         = []; % If not empty, length == 2, choose at least min and up to max based on KW test 
            defaultParams.selectFeaturesSigThresh      = 0.001; % p < this
            defaultParams.excludeFeatures              = []; % Exclude these from used features
            defaultParams.maxComponents                = []; % If empty, maxComponents will be set in .Train() based on data
            defaultParams.initProbs                    = []; % Set init probabilities
            defaultParams.pureHMM                      = 0; % If true, compute hmm through em
            defaultParams.eventInds.train              = [];
            defaultParams.eventInds.test               = [];
            defaultParams.nTimeDim                     = 0;
            defaultParams.nTimeDimStride               = 5; 
            defaultParams.xValVerbosity                = 0; % 1 = print iter, > 1 print performance
            defaultParams.xValSeed                     = []; % If not empty, seed for consistent fold inds
            defaultParams.kFold                        = 5; % Number of folds (splits) for cross validation
            defaultParams.kFoldIter                    = []; % If empty, will be set to kFold. If not empty, will iterate through kFold until we have run kFoldIter train/test xVals
            defaultParams.kFoldShuffleTrainLabels      = false; % Set true to test a "null" decoder
            defaultParams.xvalThreshold                = [];
            defaultParams.trialBoundaries              = []; % see GetTrialBounds. Used for crossval to divide folds into trials [trialStartInd trialStopInd] size [nTrial 2]
            defaultParams.labelPerTrial                = []; % Used for evenly splitting up crossval. If empty, we use the start of each trial for the label
            defaultParams.smoothDataFun                = @movmean; %'ExpSmooth'; % @movmean; % movmean, GaussSmooth2. For smoothing the input and test data
            defaultParams.smoothDataAmount             = [0 0]; % Additional param to smoothDataFun [nBack nForward]
            
            defaultParams.smoothLikFun                 = @movmean; %'ExpSmooth'; % @movmean; % movmean, GaussSmooth2. For smoothing the likelihoods
            defaultParams.smoothLikAmount              = [0 0]; % Additional param to smoothLikFun [nBack nForward]
            
            defaultParams.enforceTransitionMatrix      = []; % If empty, uses empirical trans matrix. If scalar, sets trans matrix diag to scalar value and other 1-scalar. If matrix, sets to matrix.
            defaultParams.enforceProjection            = []; % If empty, coefficients.projection uses dimReductionMethod. Otherwise uses this value
            defaultParams.contrainStateCov             = ''; % Constraints on state covariance estimates, see dHMM.m
            
            defaultParams.transform                    = []; % Empty or procrustes-like transform struct
            
            
            defaultParams.summary.projFigNum           = []; % If empty, will create a new figure each time;
            defaultParams.summary.likelihoodFigNum     = []; % If empty, will create a new figure each time;
            defaultParams.summary.confusionFigNum      = []; % If empty, will create a new figure each time;
            
            defaultParams.onSave.keepDataAndLabels     = false; % If true, we save data and labels when this object is saved (to a .mat file)
            
            if nargin < 1
                params = [];
            end
            if nargin < 2
                data = [];
            end
            if nargin < 3
                labels = [];
            end
            [params,data,labels,sSessionStruct] = CheckInputs(obj,params,data,labels);
            
            
            
            
            obj.params = MergeBintoA(defaultParams,params);
            obj.data   = data;
            obj.labels = labels;
            CheckDependencies(obj);
            
            
            % Update coefficients / params if session struct was passed
            switch sSessionStruct.passedStruct
                case ''
                    % Do nothing
                case 'slc'
                    obj.UpdateCoefficientsFromSLC(sSessionStruct.sessionStruct);
                case 'singleBlock'
                    obj.UpdateCoefficientsFromSingleBlock(sSessionStruct.sessionStruct);
                case 'sFilter'
                    obj.UpdateCoefficientsFromSFilters(sSessionStruct.sessionStruct);
                case 'dMultistate'
                    obj.UpdateCoefficientsFromMultistate(sSessionStruct.sessionStruct);
                otherwise
                    error('Unknown session object %s', sSessionStruct.passedStruct)
            end
        end
        
        % -----------------------------------------------------------------
        %   Train 
        % -----------------------------------------------------------------
        function Train(obj, trainInds)
        % Calibrate the Multistate decoder
            if nargin < 2 
                trainInds = [];
            end
            if isempty(trainInds)
                trainInds = obj.params.eventInds.train;
            end
            % If no train inds, train on all data
            if isempty(trainInds)
                trainInds = 1:size(obj.data,1);
            end
            
            if ~iscategorical(obj.labels)
                obj.labels = categorical(obj.labels);
            end
            
            if isempty(obj.params.maxComponents)
                obj.params.maxComponents = length(categories(obj.labels))-1;
            end
            
            % Find data projection
            proj = DimReduct(obj, trainInds);
            
            % Compute HMM params
            hmm        = dHMM(obj.params);
            hmm.data   = obj.info.train.projData;
            hmm.labels = obj.info.train.labels;
            hmm.Train();

            % Update coefficients, info
            obj.info.train.empiricalTransitionMatrix = hmm.info.train.empiricalTransitionMatrix;            
            obj.coefficients.projection  = proj;
            obj.coefficients             = obj.CopyBtoA(obj.coefficients,hmm.coefficients);
        end
        
        
        % -----------------------------------------------------------------
        %   Test
        % -----------------------------------------------------------------
        function [predState, liklihoods] = Test(obj, testData, testLabels)
            % Test your data
            %
            % Compute the following metrics:
            % auc
            % confustion matrix
            % precision
            % recall
            % F-measure


            % Pass in both testData and testLables,
            % Otherwise, use training data and labels

            if nargin < 2
                testInds   = obj.params.eventInds.test;
                testData   = obj.data;
                testLabels = obj.labels;
            else
                testInds = [];
            end

            if isempty(testInds)
                testInds = 1:size(testData,1); 
            end
            if ~iscategorical(testLabels)
                testLabels = categorical(testLabels);
            end

            [testData, testLabels] = PreProcess(obj, testData, testLabels, testInds);

            % Make sure test labels match train label order
            % This is done in dHMM
%             testLabels = MakeCatCommon(obj.info.train.labels, testLabels);
            

            projData = ProjectData( obj, testData );
            obj.info.test.projFeatures = projData;

            % Compute likelihoods
            hmm = dHMM(obj.params);
            hmm.coefficients = obj.coefficients;
            hmm.labels = categorical(obj.info.train.uStates, obj.info.train.uStates ); % To ensure we can make common with test labels
            hmm.Test( projData, testLabels );



            % Save off metrics and info from the hmm
            obj.info.test   = obj.CopyBtoA(obj.info.test,hmm.info.test);
            obj.performance = obj.CopyBtoA(obj.performance,hmm.performance);


            if nargout > 0
                predState  = obj.info.test.predState;
                liklihoods = obj.info.test.estLiks;
            end
        end
        
        % -----------------------------------------------------------------
        %   Predict 
        % -----------------------------------------------------------------
        function [hmmLiks, obsLiks] = Predict(obj, testData)
            % Compute liklihoods from data without labels
            if nargin < 2
                testInds = obj.params.eventInds.test;
                testData = obj.data;
            else
                testInds = [];
            end
            if isempty(testInds)
                testInds = 1:size(testData,1); 
            end
            testData = PreProcess(obj, testData, [], testInds);
            projData = ProjectData( obj, testData );

            % Compute likelihoods
            hmm = dHMM( obj.params );
            hmm.info.train.empiricalTransitionMatrix = obj.info.train.empiricalTransitionMatrix;
            hmm.coefficients = obj.coefficients;
            [hmmLiks, obsLiks] = hmm.Predict( projData );
        
        end
        
        % -----------------------------------------------------------------
        %   Cross validate
        % -----------------------------------------------------------------
                
        function xVal(obj)
            LocalxVal(obj);
        end
        
        % -----------------------------------------------------------------
        %   Get Test data likelihoods
        % -----------------------------------------------------------------
        function [liks, obsLiks] = GetMostLikely(obj,data)
            % Compute liklihoods from data without labels
            if nargin < 2
                testInds = obj.params.eventInds.test;
                testData = obj.data;
            else
                testData = data;
                testInds = [];
            end
            if isempty(testInds)
                testInds = 1:size(testData,1); 
            end
            testData = PreProcess(obj, testData, [], testInds);
            projData = ProjectData( obj, testData );

            % Compute likelihoods
            hmm = dHMM( obj.params );
            hmm.info.train.empiricalTransitionMatrix = obj.info.train.empiricalTransitionMatrix;
            hmm.coefficients = obj.coefficients;
            [liks, obsLiks] = hmm.Viterbi( projData );
        end
        
        % -----------------------------------------------------------------
        %   Get projection
        % -----------------------------------------------------------------
        function projData = GetProjection(obj, testData)
            % Compute projection based on trained dim reduct method 
            if nargin < 2
                testInds = obj.params.eventInds.test;
                testData = obj.data;
            else
                testInds = [];
            end
            testData = PreProcess(obj, testData, [], testInds);
            projData = ProjectData( obj, testData );
        end 
        
        % Function for saving decoder

        % -----------------------------------------------------------------
        %   Summarize performance 
        % -----------------------------------------------------------------
        function handles = Summarize(obj)
            % Plot/summarize/save results
            handles.proj_h = PlotCloud(obj);        
            handles.prob_h = PlotTestProbability(obj);
            handles.conf_h = PlotConfusion(obj);

            if ~nargout
                clear handles;
            end
        
        end
        
        % Function for saving decoder
        
        
    end
    
    
    %% Utils 
    
    methods 
        function RemoveDataAndLabels(obj)
        % Remove data and labels to prevent hogging memory
            obj.data = [];
            obj.labels = [];
        end
        
        function ClearDataAndLabels(obj)
        % Remove data and labels to prevent hogging memory
            RemoveDataAndLabels(obj);
        end
        
        function DeleteDataAndLabels(obj)
            RemoveDataAndLabels(obj);
        end
        
        function SelectTopFeatures(obj, trainInds)
            % Sets param featInds based on criteria (significant features)
            obj.params.featInds = 1:size(obj.data,2);
            [tmpFeat, tmpLabels] = PreProcess(obj, obj.data, obj.labels, trainInds);
            
            [selectedFeatInds, rankFeat] = SelectSigFeatures(tmpFeat, tmpLabels, ...
                    obj.params.excludeFeatures, ...
                    obj.params.selectFeaturesMinMax, ...
                    obj.params.selectFeaturesSigThresh);
            obj.params.featInds = selectedFeatInds;
            
            if isempty(selectedFeatInds)
                obj.params.selectFeaturesMinMax = [];
            end
            
            obj.info.featureSelection = rankFeat;
        end
        
        
        function UpdateCoefficientsFromSingleBlock(obj, singleBlock, useMultistate2)
            if nargin < 3
                useMultistate2 = false;
            end
            
            slc.sSLC = singleBlock.sSLCsent;
            slc.task.goalState.header.actionMap = singleBlock.sRTC.ActiveDiscreteFilter.buildDX.actionMap;            
            UpdateCoefficientsFromSLC(obj, slc, useMultistate2)
        end
        
        function UpdateCoefficientsFromSFilters(obj, sFilter, useMultistate2)
            if nargin < 3
                useMultistate2 = false;
            end
            multistateStr = 'multistate';
            if useMultistate2
                multistateStr = 'multistate2';
            end
            
            if length(sFilter) > 1
                error('Please pass a single instance of sFILTERS. E.g. sFILTERS(3) for filter 3.')
            end
            
            slc.task.goalState.header.actionMap = sFilter.buildDX.actionMap;
            
            error('Not fully supported. Need to map to sSLC')
            
            
            %% Need to pass these to sSLC
            m = sFilter.buildDX.mGesture;
            slc.sSLC.decoders.(multistateStr).enable = 1;
            
            slc.sSLC.decoders.(multistateStr).mu;
            slc.sSLC.decoders.(multistateStr).cov;
            slc.sSLC.decoders.(multistateStr).projection;
            slc.sSLC.decoders.(multistateStr).transitionMatrix;
            slc.sSLC.decoders.(multistateStr).maxNumPCs;
            
            
            sCh = GetChannelsUsed( slc );
            onlineP.feat = sCh.(multistateStr);
            featFN = fieldnames( onlineP.feat.ch );
            smoothSubfield = [featFN{1} 'SmoothSteps'];
            slc.sSLC.decoders.(multistateStr).(smoothSubfield) = m.params.smoothDataAmount
            slc.sSLC.decoders.postProcess.discrete.smoothSteps;

            
            
            UpdateCoefficientsFromSLC(obj, slc, useMultistate2)
        end
        function UpdateCoefficientsFromMultistate(obj,m)
            props = properties(m);
            if isempty(props)
                props = fieldnames(m); % May have been converted to struct
            end
            props(ismember(props, 'params')) = []; % Already passed, merged
            for ii = 1:length(props)
                obj.(props{ii}) = m.(props{ii});
            end
            
        end
        
        function UpdateCoefficientsFromOnline(obj, slc, useMultistate2)
            newFunName = 'UpdateCoefficientsFromSLC';
            errorStr = sprintf('\nFunction name change. Please update code to call new function: \n%s\n', newFunName);
            error(errorStr);
        end
        
        
        function UpdateCoefficientsFromSLC(obj, slc, useMultistate2)
            % slc should have .sSLC and .task (for uStates)
            if nargin < 2
                slc = [];
            elseif isfield( slc, 'sSLCversion' )
                slc.sSLC = slc;
            end
            if nargin < 3
                useMultistate2 = false;
            end
            multistateStr = 'multistate';
            if useMultistate2
                multistateStr = 'multistate2';
            end
            
            

            if isempty(slc)
                % Return to original params / coefficients
                if isfield( obj.info, 'org')        
                    obj.params.enforceTransitionMatrix = [];
                    obj.coefficients = obj.info.org.coefficients;
                    obj.params = MergeBintoA( obj.params, obj.info.org.params );
                    obj.info.train.empiricalTransitionMatrix = obj.info.org.info.train.empiricalTransitionMatrix;
                    obj.info = rmfield(obj.info, 'org');
                end
            else
                % Update with SLC
                if ~isfield( obj.info, 'org')
                    obj.info.org.coefficients = obj.coefficients;
                    obj.info.org.params = obj.params;
                    % Maybe save off all info?
                    try % Try to set, otherwise, set to empty?
                        obj.info.org.info.train.empiricalTransitionMatrix = obj.info.train.empiricalTransitionMatrix;
                    catch
                        obj.info.org.info.train.empiricalTransitionMatrix = [];
                    end
                end
                
                if isfield(slc,'task')
                    uStates = unique({slc.task.goalState.header.actionMap.imagery}', 'stable');
                else
                    uStates = {};
                end


                if isfield( slc.sSLC.decoders, 'PCAHMM') && slc.sSLC.decoders.PCAHMM.enable
                    sCh = GetChannelsUsed( slc );

                    onlineP.feat = sCh.PCAHMM;
                    %%
                    featFN = fieldnames( onlineP.feat.ch );
                    smoothSubfield = [featFN{1} 'SmoothSteps'];
                    onlineP.smooth = slc.sSLC.decoders.PCAHMM.(smoothSubfield);
                    onlineP.smoothLik = slc.sSLC.decoders.postProcess.discrete.smoothSteps;
                    %%
                    onlineP.emissionMu = slc.sSLC.decoders.PCAHMM.mu;
                    onlineP.emissionCov = slc.sSLC.decoders.PCAHMM.cov;
                    onlineP.projection = slc.sSLC.decoders.PCAHMM.projection;
                    onlineP.transitionMatrix = slc.sSLC.decoders.PCAHMM.transitionMatrix;
                    onlineP.numComponents = slc.sSLC.decoders.PCAHMM.maxNumPCs;
                    onlineP.numOutStates = slc.sSLC.decoders.PCAHMM.maxStates;
                    

                elseif isfield( slc.sSLC.decoders, multistateStr) && slc.sSLC.decoders.(multistateStr).enable
                    sCh = GetChannelsUsed( slc );
                    onlineP.feat = sCh.(multistateStr);
                    featFN = fieldnames( onlineP.feat.ch );
                    smoothSubfield = [featFN{1} 'SmoothSteps'];
                    onlineP.smooth = slc.sSLC.decoders.(multistateStr).(smoothSubfield);
                    onlineP.smoothLik = slc.sSLC.decoders.postProcess.discrete.smoothSteps;

                    onlineP.emissionMu = slc.sSLC.decoders.(multistateStr).mu;
                    onlineP.emissionCov = slc.sSLC.decoders.(multistateStr).cov;
                    onlineP.projection = slc.sSLC.decoders.(multistateStr).projection;
                    onlineP.transitionMatrix = slc.sSLC.decoders.(multistateStr).transitionMatrix;
                    onlineP.numComponents = slc.sSLC.decoders.(multistateStr).maxNumPCs;
                    onlineP.numOutStates = slc.sSLC.decoders.(multistateStr).maxStates;
                else
                    error('Did not find PCAHMM or Multistate decoder. Was a discrete decoder enabled?')
                end

                sSlcFeatInds = 1:length(onlineP.feat.featInds);
                mxS = onlineP.numOutStates;
                if isempty(obj.params.maxComponents)
                    mxC = onlineP.numComponents;        
                    obj.params.maxComponents = mxC;
                else
                    mxC = min( obj.params.maxComponents, onlineP.numComponents );
                end
                obj.info.train.empiricalTransitionMatrix(1:mxS,1:mxS) = onlineP.transitionMatrix(1:mxS,1:mxS);
                obj.info.train.uStates = uStates;
                obj.params.enforceTransitionMatrix(1:mxS,1:mxS) = onlineP.transitionMatrix(1:mxS,1:mxS);
                obj.params.smoothDataAmount = [onlineP.smooth 0];
                obj.params.smoothLikAmount = [onlineP.smoothLik 0];
                obj.params.featInds = onlineP.feat.featInds;
                obj.coefficients.transitionMatrix(1:mxS,1:mxS) = onlineP.transitionMatrix(1:mxS,1:mxS);
                obj.coefficients.projection = onlineP.projection(sSlcFeatInds,1:mxC);
                obj.coefficients.emissionMu(1:mxS,1:mxC) = onlineP.emissionMu(1:mxS,1:mxC);
                obj.coefficients.emissionCov(1:mxS,1:mxC, 1:mxC) = onlineP.emissionCov(1:mxS,1:mxC, 1:mxC);
                

            end
        end
        
        %% Cross validation
        % -----------------------------------------------------------------
        %   Cross validate
        % -----------------------------------------------------------------
        function LocalxVal(obj)
        
            if ~isempty(obj.params.kFoldIter)
                kFoldIter = obj.params.kFoldIter;
            else
                kFoldIter = obj.params.kFold;
            end
            if ~isempty(obj.params.xValSeed)
                orgRng = rng();
                rng(obj.params.xValSeed);
            end
            
            metric  = nan(kFoldIter,1);
            
            % Remove previous xvalidation entries
            if isstruct(obj.performance)
                if isfield(obj.performance, 'xPerformance')
                    obj.performance = rmfield(obj.performance, 'xPerformance');
                end
                if isfield(obj.performance, 'std')
                    obj.performance = rmfield(obj.performance, 'std');
                end
            end
            
            if obj.params.kFoldShuffleTrainLabels
            end
            
            
            %% Perform cross validation
            iter = 1;
            done = 0;
            startStops = GetStartStopsFromTrialBounds(obj);
            while ~done
                
                cvInds = GetFoldInds(obj);
                
                for kk = 1:obj.params.kFold
                    if obj.params.xValVerbosity
                        fprintf('XVal %d of %d\n', iter, kFoldIter);
                    end
                    
                    testInds  = cvInds == kk;
                    trainInds = find(~testInds);
                    
                    obj.Train(trainInds);

                    [~, liks] = obj.Test( obj.data(testInds,:), obj.labels(testInds,:) );
                    %%
                    if ~isempty(startStops)
                    testLabels = obj.labels(testInds,:);
                    testStStp = startStops(cvInds(startStops(:,1))==kk,:);
%                     threshold = max(obj.performance.rocThresh);
                    
                    threshold = 0.95;
                    if ~isempty(obj.params.xvalThreshold)
                        threshold = obj.params.xvalThreshold;
                    end
                    if isnumeric(threshold)
                        prctFirstCorrect = ComputePercentFirstCorrect(obj,liks,testLabels,testStStp, threshold);
                    elseif strcmp(threshold,'best')
                        threshs = [0.85 0.9 0.95 0.97 0.99 0.999 0.9999 0.99999];
                        for ii = 1:length(threshs)
                            threshold = threshs(ii);
                            pfc(ii) = ComputePercentFirstCorrect(obj,liks,testLabels,testStStp, threshold);
                        end
                        [mv, mi] = max(pfc);
                        fprintf('Best threshold %.6f with %.3f\n', threshs(mi), mv);
                        prctFirstCorrect = mv;
                    end
                    else
                        prctFirstCorrect = nan;
                    end
                    %%

                    metric(iter) = obj.performance.auc(2);
                    perf = obj.performance;
                    perf.prctFirstCorrect = prctFirstCorrect;
                    xPerformance(iter) = perf;
                    
                    if obj.params.xValVerbosity > 1
                        fprintf('   Perf: avg AUC: %.3f  MCC: %.3f\n', mean(obj.performance.auc(2:end)), obj.performance.mcc);
                    end
                    
                    iter = iter + 1;
                    done = iter > kFoldIter;
                    if done
                        break
                    end

                end
            end
            
            AverageXValPerformance(obj, xPerformance);
            
            % Restore the random number generator
            if ~isempty(obj.params.xValSeed)
                rng(orgRng);
            end

        end
        
        function AverageXValPerformance(obj, xPerformance)
            % Average the performance
            % For each performance metric, concatenate the metric then
            % average.
            fn = fieldnames(xPerformance);
            fn(ismember(fn, {'xPerformance', 'std'})) = []; % remove xPerformance if it exists
            for ii = 1:length(fn)
                % First get size of dims to cat on first singleton dim
                sz = size( xPerformance(1).(fn{ii}) );
                catInd = find(sz==1,1);
                if isempty(catInd)
                    catInd = length(sz)+1; % If no singleton dim, concat on the next dim
                end
                catPerf = cat( catInd, xPerformance.(fn{ii}) );
                
                % Average over the cross validated results
%                 noNanPerf = catPerf(~isnan(catPerf));
                obj.performance.(fn{ii}) = nanmean(catPerf, catInd);
                obj.performance.std.(fn{ii}) = nanstd(catPerf, [], catInd);
            end
            
            % Save the raw cross validated performance
            obj.performance.xPerformance = xPerformance;
        end
        
        function [prctFirstCorrect,isFirstCorrect, firstPredInds] = ComputePercentFirstCorrect(obj,liks,testLabels,testStStp, threshold)
            
            
            [maxV,maxI] = max(liks,[],2);
            if ~isempty(threshold)
                maxI(maxV<threshold) = 1;
            end

            trialBounds = GetTrialBounds(testStStp, 'startStops', 1);
            nTestTrials = length(trialBounds)-1;
            isFirstCorrect = nan(nTestTrials,1);
            firstPredInds= nan(nTestTrials,1);
            for ii = 1:nTestTrials
                testTrlInds = trialBounds(ii):trialBounds(ii+1)-1;
                testTrlLbls = double(testLabels(testTrlInds));
                if any(testTrlLbls==1) %any(ismember(testTrlLbls,'no_action'))
                    continue
                end
                testTrlPred = maxI(testTrlInds);
                firstPredInd = find(testTrlPred ~= 1,1);
                if ~isempty(firstPredInd)
                    isFirstCorrect(ii) = testTrlPred(firstPredInd) == testTrlLbls(1);
                    firstPredInds(ii) = firstPredInd;
                else
                    isFirstCorrect(ii) = 0;
                end
            end

            prctFirstCorrect = sum(isFirstCorrect == 1)/sum(~isnan(isFirstCorrect))*100;
            if obj.params.xValVerbosity > 1
                fprintf('%.3f first correct\n', prctFirstCorrect)
            end

        end
        function startStops = GetStartStopsFromTrialBounds(obj,trialBoundaries)
            if nargin < 2 || isempty(trialBoundaries)
                trialBoundaries = obj.params.trialBoundaries;
            end
                
            if isvector(trialBoundaries)
                % assume trialBoundaries are the start inds
                % final index is the equivalent to the stop trial index
                % + 1
                for ii = 1:length(trialBoundaries)-1
                    startStops(ii,:) = [trialBoundaries(ii), trialBoundaries(ii+1)-1];
                end

            else
                startStops = trialBoundaries;
            end
        end
        
        
        function cvInds = GetFoldInds(obj)
        % Gets fold indices
        % Expects the following fields / parameters
        % Required:
        %   obj.params.kFold
        %   obj.data
        %
        % Optional: (if ~isempty(obj.params.trialBoundaries))
        %   obj.params.trialBoundaries - shuffle trials instead of points
        %   obj.params.eventInds.train
            nPoints = size(obj.data,1);
            
            
            
            if ~isempty(obj.params.trialBoundaries)
                trialStartStops = GetStartStopsFromTrialBounds(obj);
                
                
                % Select trials that coincide with the eventInds
                if ~isempty( obj.params.eventInds.train )
                    validTrials = ismember( trialStartStops, obj.params.eventInds.train );
                    trialStartStops = trialStartStops(validTrials(:,1), :);
                end
                
                
                nTrials = size(trialStartStops,1);
                kFold = obj.params.kFold;
                % Divide the trials into the folds

                %%
                try
                    if ~isempty(obj.params.labelPerTrial) && length(obj.params.labelPerTrial) == size(trialStartStops,1)
                        % User specified labels per trial
                        lblPerTrial = obj.params.labelPerTrial;
                    else
                        lblPerTrial = obj.labels(trialStartStops(:,1),:);
                    end
                    
                catch
                    error('Cross validation error when indexing trials into labels.\n')
                    %keyboard
                end
                uLbls = unique(lblPerTrial);
                cvTrls = nan(nTrials,1);
                for ii = 1:length(uLbls)
                    lblTrlsInds = find(lblPerTrial == uLbls(ii));
                    nTrialInds = length(lblTrlsInds);
                    nPerFold = floor( nTrialInds/kFold );
                    if nPerFold == 0
                        % Not a full set, set the first N folds
                        foldNumber = 1:nTrialInds;
                        cvTrls(lblTrlsInds) = foldNumber;
                    else
                        for kk = 1:kFold
                            foldInds = lblTrlsInds( randperm(length(lblTrlsInds), nPerFold) ); 
                            lblTrlsInds = setdiff(lblTrlsInds,foldInds);
                            cvTrls(foldInds) = kk;
                        end
                        % Any remaining inds get distributed randomly
                        cvTrls(lblTrlsInds) = randperm(kFold,length(lblTrlsInds));
                    end
                    
                end
                %%
%                 cvTrls = crossvalind('Kfold',nTrials,kFold);
                
                %% Set each trial to it's associated fold number
                cvInds = nan(nPoints,1);
                for kk = 1:kFold
                    selTrials = ismember( cvTrls, kk );
                    trialInds = RowColon( trialStartStops(selTrials,:) );
                    cvInds(trialInds) = kk;
                end
            
            else
                cvInds = crossvalind('Kfold',nPoints,obj.params.kFold);
            end
        end
    end
    
    
    
    %% Plot
    methods
        function handles = PlotCloud(obj)
%             projData = obj.data*obj.coefficients.projection;
            pltData = {};
            pltLbls = {};
            pltTitles = {};
            if ~isempty( obj.info )
                
                
                if isfield(obj.info, 'train')
                    pltData{end+1}   = obj.info.train.projData;
                    pltLbls{end+1}   = double(obj.info.train.labels);
                    pltTitles{end+1} = 'Train projection';
                    uStates = obj.info.train.uStates;
                end
                
                if isfield(obj.info, 'test')
                    pltData{end+1}   = obj.info.test.projFeatures;
                    pltLbls{end+1}   = double(obj.info.test.labels);
                    pltTitles{end+1} = 'Test projection';
                    uStates = obj.info.train.uStates;
                end
            end
            
            if isempty(pltData)
                if nargout
                    handles = gobjects(0);
                end
                return;
            end
            
            
            nPlts   = length(pltData);
            nDim    = size(pltData{1},2);
            clrs    = obj.GetColors( length(uStates) );            
            clrs(end+1,:) = ones(1,3).*0.3;
            nanIndex = size(clrs,1);
            
            
            
            if ~isempty(obj.params.summary.projFigNum)
                fig_h = figure( obj.params.summary.projFigNum );
            else
                fig_h = figure;
            end
            clf;
            
            anyNans = false;
            for spI = 1:nPlts
                if nPlts > 1
                    ax(spI) = subplot(1,nPlts,spI);
                else
                    ax = gca;
                end
                
                % Update nan labels
                nanLogical = isnan(pltLbls{spI});                
                pltLbls{spI}(nanLogical) = nanIndex;
                anyNans = anyNans || any(nanLogical);
                
                
                if nDim == 1
                    for ii = 1:length(uStates)
                        selI = pltLbls{spI} == ii;
                        [fout0,xout] = ksdensity( pltData{spI}(selI) );
                        plot(xout,fout0, 'color', clrs(ii,:));
                    end
                elseif nDim == 2
                    scatter(ax(spI), pltData{spI}(:,1),  pltData{spI}(:,2), 5, clrs(pltLbls{spI},:), 'filled')
                else
                    
                    scatter3(ax(spI), pltData{spI}(:,1),  pltData{spI}(:,2), pltData{spI}(:,3), 5, clrs(pltLbls{spI},:), 'filled');
                    view(ax(spI),3)
                end
                
                title(ax(spI), pltTitles{spI})
            end
            
            if nDim == 1
                linkaxes(ax,'xy')
            end
            
            % Plot for legend
            if anyNans
                legStr = [uStates(:); {'NaN label'}];                
            else
                legStr = uStates;
            end
            
            hold on
            for ii = 1:length(legStr)
                lph(ii) = plot(nan, 'o', 'color', clrs(ii,:), 'MarkerFaceColor', clrs(ii,:));
            end
            legend(lph, legStr, 'Location', 'best')
            
            handles.fig_h = fig_h;
            handles.axes  = ax;
            
        
        end
        
        
        function handles = PlotTestProbability(obj)
            
            if ~isfield( obj.info, 'test' )
                handles = [];
                return
            end
            
            grpInds = double( obj.info.test.labels );
            uStates = obj.info.train.uStates;
            estLiks = obj.info.test.estLiks;
            predState = double(obj.info.test.predState);
            incorrectGuess = obj.info.test.incorrectGuess;
            nTrainedStates = length(uStates);
            nPoints = length(grpInds);
            nDecoderStates = size(estLiks,2);
            
            
            annPos = [0.01 0.9, 0.08 0.05];
            annTxtYOffset = 0.02;
            yOffset = 0.1;
            yLabelOffset = 1 + yOffset/3;
            yGuessOffset = 1 + yOffset/2;
            yIncrctOFfset = 1 + yOffset * 2/3;
            xLabelOffset = -round(nPoints*0.1);
            
            if ~isempty(obj.params.summary.likelihoodFigNum)
%                 alreadyFig = ~isempty(findobj('Type', 'figure', 'number', obj.params.summary.likelihoodFigNum));
                fig_h = figure( obj.params.summary.likelihoodFigNum );
            else
%                 alreadyFig = false;
                fig_h = figure;
            end
            
            fig_h.Color = [1 1 1];
            clf
            
            
            %--------------------------------------------------------------
            % Probability axes
            %--------------------------------------------------------------
            axMain = gca;
            hold on;
            % plot(labels, 'k', 'LineWidth', 3)
            clrs = obj.GetColors( nTrainedStates );
            
            incrt_h = gobjects(0);
            for ii = 1:nTrainedStates
                if ii <= nDecoderStates
                    plot(estLiks(:,ii), '-', 'color', clrs(ii,:), 'LineWidth', 0.3);
                end
                bgInd = find(predState==ii);
                pltVals = nan(nPoints,1);
                pltVals(bgInd) = yGuessOffset;
                plot(pltVals, 'Color', clrs(ii,:), 'LineWidth', 3)

                pltVals = nan(nPoints,1);

                pltVals(grpInds==ii) = yLabelOffset;
                ph(ii) = plot(pltVals, 'Color', clrs(ii,:), 'LineWidth', 5);

                pltIncrtInd = find( incorrectGuess(:) & predState==ii );
                if ~isempty(pltIncrtInd)
                    incrt_h = plot(pltIncrtInd, ones(size(pltIncrtInd)).*yIncrctOFfset, 's', 'Color', 'r', 'markerfacecolor','r', 'MarkerSize', 1);
                end

            end
            predLblStrOpt = {'HorizontalAlignment', 'center', 'LineStyle', 'none', 'FontWeight', 'bold'};
            annotation('textbox', annPos, 'String','Incorrect', predLblStrOpt{:})
            annPos(2) = annPos(2) - annTxtYOffset;
            annotation('textbox', annPos, 'String','Pred', predLblStrOpt{:})
            annPos(2) = annPos(2) - annTxtYOffset;
            annotation('textbox', annPos, 'String','Label', predLblStrOpt{:})
            % text(xLabelOffset,yIncrctOFfset, 'Incorrect', 'HorizontalAlignment', 'center')
            % text(xLabelOffset,yGuessOffset, 'Pred', 'HorizontalAlignment', 'center')
            % text(xLabelOffset,yLabelOffset, 'Label', 'HorizontalAlignment', 'center')
            axMain.XLim(2) = length(pltVals);
            axMain.YLim(2) = yIncrctOFfset;
            axMain.YTick(axMain.YTick>1) = [];
            drawnow
            if isnumeric( uStates )
                legendStr = sprintfc('State %d', uStates);
                xTickRotation = 0;
            else
                legendStr = uStates;
                xTickRotation = 35;
            end
            legendStr{end+1} = 'Incorrect';
            [~, objh] = legend([ph incrt_h], legendStr{:}, 'location','northeastoutside', 'Interpreter', 'none');
            lineh = findobj(objh,'type','line');
            set(lineh, 'linewidth', 2);
            
            xlabel('Analysis points')
            ylabel('Probability')
            title('Multistate probability and best state guess')
            
            handles.fig = fig_h;
            handles.ax = axMain;
        end
        
        
        function handles = PlotConfusion(obj)
            
            if ~isfield( obj.info, 'test' )
                handles = [];
                return
            end
            
            uStates = obj.info.train.uStates;
            nTrainedStates = length(uStates);

            %--------------------------------------------------------------
            % Confusion Matrix
            %--------------------------------------------------------------    
            if ~isempty(obj.params.summary.confusionFigNum)
                fig_conf = figure( obj.params.summary.confusionFigNum );
            else
                fig_conf = figure;
            end
            fig_conf.Color = [1 1 1];
            
            confVal = round(obj.performance.conf(1:nTrainedStates,1:nTrainedStates)); % Round because cross validation may give non-integer values
            
            cm = confusionchart(fig_conf, confVal,categorical(uStates,uStates));
            cm.RowSummary = 'row-normalized';
            cm.ColumnSummary = 'column-normalized';
            
            prctCorrect = confVal./sum(confVal,2);
            prctCorrect = 100*mean(diag(prctCorrect));
            title(sprintf('Average correct: %.2f%%', prctCorrect))
            

            handles.fig = fig_conf;
            handles.ax = cm;
            
        end
        
    end
    
    methods 
        function timeDimData = AddTimeAsFeatDim( obj, d )
            win         = -obj.params.nTimeDim:obj.params.nTimeDimStride:0;
            dEpoch      = GetWindowAroundEvents(1:size(d,1), win)';
            nNewFeat    = size(d,2)*size(dEpoch,2);
            timeDimData = reshape( d(dEpoch,:), size(dEpoch,1), nNewFeat );
            % Pad the first n time steps with 0s (should be removed downstream)
            timeDimData = cat(1, zeros(obj.params.nTimeDim, nNewFeat), timeDimData);
        end
        
        function inds = UpdateEpochIndsFromTimeAsFeatDim(obj, inds)
            if obj.params.nTimeDim > 0
                % Removes the first nTimeDim points from each consecutive
                % epoch in inds. This will prevent the data, labels
                % from being contaminated with data outside of inds.
                % (though, currently smoothing may contaminate the data...)
                % Assumes inds are in order (monotonic)
                indDiff = diff(inds(:));
                if any(indDiff < 1)
                    error('Expected inds to be monotonic')
                end
                % Find the starts of each epoch in inds                
                indBreaks = [find(indDiff>1)+1];
                % Remove the first nTimeDim from each start
                % rmvEpoch is NSteps x nBreaks
                rmvEpoch = GetWindowAroundEvents(indBreaks,0:obj.params.nTimeDim-1);
                % Remove inds
                inds(rmvEpoch(:)) = [];
            end
        end
        
        function [data, labels] = PreProcess(obj, data, labels, inds)
            % Assumes data is already zscored
            % Adds time as a feature dimension
            % Smooths data
            %
            % Updates inds (and assumes that no other part of the code
            % needs to know about the updated inds)
            
            if ischar(obj.params.smoothDataFun)
                if ismember( obj.params.smoothDataFun, methods(obj) )
                    % Assume we are calling the local smooth function
                    smoothDataFun = @(varargin) obj.(obj.params.smoothDataFun)(varargin{:});
                else
                    % Assume it is a regular function
                    smoothDataFun = str2func( obj.params.smoothDataFun );
                end
            else
                smoothDataFun = obj.params.smoothDataFun;
            end
            
            featInds = obj.params.featInds;
            if isempty(featInds)
                featInds = 1:size(data,2);
            end
            
            featInds = setdiff(featInds, obj.params.excludeFeatures);

            % This order is important! (but maybe can be improved on...)
            % 1) Select features
            % 2) Smooth on all data
            % 3) Add prev time steps as feat
            % 4) Select time inds
            data = data(:,featInds);
            data = smoothDataFun( data, obj.params.smoothDataAmount ); % Should we smooth before adding time dim. Maybe so since we kind of "lose" nTimDim time steps.
            data = AddTimeAsFeatDim( obj, data );
            inds = UpdateEpochIndsFromTimeAsFeatDim(obj, inds);
            data = data(inds,:);
            
            if ~isempty(labels)
                labels = labels(inds);
            end
            
        end
        
        
        
        function [proj, projData] = DimReduct(obj, trainInds)
            
                          
            if isempty( obj.params.featInds ) && ~isempty( obj.params.selectFeaturesMinMax )
                SelectTopFeatures(obj, trainInds);
            end
            
            [data, trainLabels] = PreProcess(obj, obj.data, obj.labels, trainInds);
            
            switch obj.params.dimReductionMethod
                case 'lda'
                    reg = obj.params.ldaRegularization;
%                     rmvInds = any(isnan(data),2);
                    
                    [proj, projData] = obj.LDA(data, trainLabels, reg );
                case 'pca'
                    %%
%                     dOrg = data;
%                     data = dOrg;
                    minMaxThresh = prctile(data, [1,95]);
                    
%                     figure; 
%                     plot(minMaxThresh')
%                     plot(data(1890,:))
                    
                    isThresh = data < minMaxThresh(1,:);
                    fillVal = isThresh.*minMaxThresh(1,:);
                    data(isThresh) = fillVal(isThresh);
                    
                    isThresh = data > minMaxThresh(2,:);
                    fillVal = isThresh.*minMaxThresh(2,:);
                    data(isThresh) = fillVal(isThresh);
                    
                    
                    %%
                    [proj,projData]  = pca(data, 'NumComponents', obj.params.maxComponents);
                case 'dca'
                    if size(obj.data,2) ~= 384
                        error('DCA currently is expecting 2 features. Need to update to provide information on feature boundaries')
                    end
                    % Use cca to align each feature, then perform lda
                    reg = obj.params.ldaRegularization;
                    f1 = 1:192;
                    f2 = 193:size(obj.data,2);
                    [A,B,r,U,V] = canoncorr(data(:,f1),data(:,f2));
                    obj.coefficients.projA = A;
                    obj.coefficients.projB = B;
                    nd = cat(2,U,V);
                    [proj, projData] = obj.LDA( nd, trainLabels, reg );
                
                    %% Or use dcaFuse
%                     [Ax,Ay,Xs,Ys] = dcaFuse(obj.data(trainInds,f1)',obj.data(trainInds,f2)',obj.labels(trainInds)');
%                     projData = cat(1, Xs,Ys)';
%                     projData = (Xs+Ys)';
%                     proj = cat(2, Ax, Ay)';
                case 'none'
                    proj = eye(size(data,2));
                    projData = data;
                otherwise
                try
                    %%
                    X = data';
%                     X = GaussSmooth2(data,150)';
%                     X = RemoveOutliers(data)';
%                     maxOut = 15;
%                     X(abs(X)>maxOut) = maxOut;
                    r = obj.params.maxComponents;
                    Q_0 = project_stiefel(randn(size(X,1),r));
    
%                     optims = {'heuristic','grassmann_trust', 'stiefel', 'grassmann_mosd', 'stiefel_trust'};

                    [ proj , fQ , info ] = run_lda( X , trainLabels , r , obj.params.dimReductionMethod , Q_0 );
                    projData = X'*proj;
                catch
                    error('Unrecognized dim reduction method %s', obj.params.dimReductionMethod)
                end
            end
    

            pcInds = 1:obj.params.maxComponents;
            nProjDim = size(proj,2);
            if nProjDim < obj.params.maxComponents
                % Pad proj and projData
                proj = cat(2, proj, zeros(size(proj,1), obj.params.maxComponents-nProjDim));
                projData = cat(2, projData, zeros(size(projData,1), obj.params.maxComponents-nProjDim));
            else
                proj = proj(:,pcInds);
                projData = projData(:,pcInds);
            end
            if ~isempty(obj.params.enforceProjection)
                projData = data*obj.params.enforceProjection;
            end
            
            if ~isempty(obj.params.transform)
                procT = obj.params.transform;
                u = mean(projData);
                projData = projData - u;
                projData = procT.b * projData * procT.T + procT.c(1,:);
                projData = projData + u;
            end
            
            % Save
            obj.info.train.uStates  = categories(trainLabels);            
            obj.info.train.labels   = trainLabels;
            obj.info.train.projData = projData;
                
        end
        
        function [eMu, eCov, transitionMatrix] = CalculateStateMeansVarsAndTransitionMatrix(obj)
            
            uStates  = obj.info.train.uStates;
            grpInd   = double(obj.info.train.labels);
            projData = obj.info.train.projData;
            

            nStates = length(uStates);

            
            transitionMatrix = nan(nStates);
            for ii = 1:nStates
                for jj = 1:nStates
                    % Select the tranistion inds from state ii to jj
                    currentStateInds = grpInd(1:end-1) == ii;
                    nextStateInds = grpInd(2:end) == jj;

                    selI = currentStateInds & nextStateInds;

                    % Count the prob of this transition for this state
                    totalTrans   = sum(selI);
                    totalInState = sum( grpInd == ii );

                    % If this state doesn't exist, set trans to 0
                    if totalInState == 0
                        transitionMatrix(ii,jj) = 0;
                    else
                        transitionMatrix(ii,jj) = totalTrans/totalInState;
                    end

                end
            end


            %% Compute emmission means, cov for each state
            nComp = obj.params.maxComponents;
            eMu   = nan(nStates, nComp);
            eCov  = nan(nStates, nComp, nComp);
            
            for ii = 1:nStates
                eMu(ii,:)    = mean(projData(grpInd == ii,:));
                eCov(ii,:,:) = cov(projData(grpInd == ii,:));
            end


        end
       
        
        function projData = ProjectData( obj, testData )
            if strcmpi(obj.params.dimReductionMethod, 'dca')
                meanData = zeros(1,384); %mean(testData);
                f1 = 1:192;
                f2 = 193:384;
                U = (testData(:,f1)-meanData(f1))*obj.coefficients.projA;
                V = (testData(:,f2)-meanData(f2))*obj.coefficients.projB;
                testData = cat(2,U,V);
            end
            % Project
            if ~isempty(obj.params.enforceProjection)
                projW = obj.params.enforceProjection;
            else
                projW = obj.coefficients.projection;
            end
            if ~isempty(projW)
                projData = testData*projW;
            else
                projData = testData;
            end
            
            if ~isempty(obj.params.transform)
                procT = obj.params.transform;
                u = mean(projData);
                projData = projData - u;
                projData = procT.b * projData * procT.T + procT.c(1,:);
                projData = projData + u;
            end
        end
        
        function probs = ComputePred( obj, testData, testLabels )
            
            
            
            
            projData = ProjectData( obj, testData );
            obj.info.test.projFeatures = projData;
            
            % Compute likelihoods
            hmm = dHMM(obj.params);
            hmm.coefficients = obj.coefficients;
            hmm.info.train.empiricalTransitionMatrix = obj.info.train.empiricalTransitionMatrix;
            
            probs = hmm.Test( projData, testLabels );


            % Save off metrics and info from the hmm
            obj.info.test = CopyBtoA(obj.info.test,hmm.info);
            obj.performance = CopyBtoA(obj.performance,hmm.performance);

        end
        

        function [params,data,labels, sSessionStruct] = CheckInputs(obj,params,data,labels)
            % If labels is empty, but params and data are not. Assume
            % data,labels were passed.
            wasParamsSkipped = ~isempty(params) && isnumeric(params) && ~isempty(data) && isempty(labels);
            if wasParamsSkipped
                labels = data;
                data = params;
                params = [];
            end
            
            if iscategorical(labels)
                % If categorical label, reset to only include categories in
                % the training data.
                % removecats removes unused categories while maintaining
                % original category order.
                labels = removecats(labels);
            end
            
            % Check if a session struct (slc, sSLC, SingleBlock, or sFILTER) was passed
            
            wasSLCPassed = isstruct(params) && isfield(params, 'sSLC');
            wasSSLCPassed = isstruct(params) && isfield(params, 'sSLCversion');
            wasSingleBlockPassed = isstruct(params) && isfield(params, 'sSLCsent') && isfield(params, 'sRTC');
            wasSFilterPassed = isstruct(params) && isfield(params, 'buildDX');
            wasMultistate = isa(params, 'dMultistate') || (isstruct(params) && isfield(params, 'coefficients'));
            
            
            if wasSLCPassed
                sSessionStruct.sessionStruct = params;
                sSessionStruct.passedStruct = 'slc';
                params = [];
            elseif wasSSLCPassed
                sSessionStruct.sessionStruct.sSLC = params;
                sSessionStruct.passedStruct = 'slc'; % slc uses sSLC
                params = [];
            elseif wasSingleBlockPassed
                sSessionStruct.sessionStruct = params;
                sSessionStruct.passedStruct = 'singleBlock';
                params = [];
            elseif wasSFilterPassed
                sSessionStruct.sessionStruct = params;
                sSessionStruct.passedStruct = 'sFilter';
                params = [];
            elseif wasMultistate
                sSessionStruct.sessionStruct = params;
                sSessionStruct.passedStruct = 'dMultistate';
                params = params.params;
                
            else
                sSessionStruct.sessionStruct = struct([]);
                sSessionStruct.passedStruct = '';
            end
        end
        
        function CheckDependencies(obj)
            % Make sure we have the required libs
            
            %% Check if HMM lib is on path
            if isempty(which('mixgauss_init'))
                hmmLibNotFound = true; % Apologies for the negative logic
                
%                 % Attempt to locate it if we are in the analysis repo
%                 classFileLoc         =  mfilename('fullpath');
%                 userAnalysisFolder   = 'UserAnalysis';
%                 isUserAnalysisOnPath = regexp(classFileLoc, userAnalysisFolder);
% 
%                 if ~isempty(isUserAnalysisOnPath)
%                     repoRoot = classFileLoc( 1:isUserAnalysisOnPath-1);
%                     expectedHMMLibPath = fullfile(repoRoot,'Shareware', 'HMM_subset');
%                     if isfolder(expectedHMMLibPath)
%                         hmmLibNotFound = false;
%                         addpath(genpath(expectedHMMLibPath))
%                     end
%                 end
                
                if hmmLibNotFound
                    error('HMM Toolbox not on path! This decoder requires Kevin Murphy''s HMM Toolbox.')
                end
                
            end
        end
        
        function s = saveobj(obj)
            if ~obj.params.onSave.keepDataAndLabels
                obj.data = [];
                obj.labels = [];
            end
            s = obj;
        end
    end
    
    methods (Static)
        function data = ExpSmooth(data,alpha)
            for ii = 2:length(data)
                data(ii,:) = alpha .* data(ii-1,:) + (1-alpha).*data(ii,:);
            end
        end
        
        function a = CopyBtoA(a,b)
            fn = fieldnames(b);
            for ii = 1:length(fn)
                a.(fn{ii}) = b.(fn{ii});
            end
        end
        
       
        
        
        function [W, Y, lambda] = LDA(X, label, reg)
        % Linear discriminate analysis with regularization
        % Inputs:
        %   X
        %       Data matrix [nPoints x nFeat]
        %
        %   label
        %       Category label vector with a length of nPoints
        %
        %   reg
        %       Regularizer parameters (scalar). Should be between 0 and 1.
        %       Shrinks the within-class scatter.
        %       See Friedman, J. H. (1989). "Regularized Discriminant Analysis"
        %
        % Outputs:
        %
        %   W
        %       Sorted eigenvectors from the solution of the LDA problem.
        %       Note, W will have a rank of nClasses-1
        %
        %   Y
        %       X projected onto W, i.e. Y = X*W;
        %
        %   lambda
        %       Eigenvalues of W
        %
        
            if nargin < 3
                reg  = 0; % 0 is no regularization.
            end
%%
            nFeat    = size(X, 2);
            [classes,~,uGrps] = grp2idx(label(:));
            nClasses = numel(uGrps);

            fullMean = nanmean(X);
            Sw       = zeros(nFeat);
            Sb       = zeros(nFeat);
            for j=1:nClasses
                Xj = X(j == classes,:);

                classMean         = nanmean(Xj);
                classFullMeanDiff = classMean-fullMean;
                classMeanCov      = (classFullMeanDiff'*classFullMeanDiff)/(size(Xj,1)-1);
                Xj                = Xj-classMean;
                withinClassCov    = (Xj'*Xj)/(size(Xj,1)-1);

                Sw = Sw + withinClassCov; 
                Sb = Sb + classMeanCov; 
            end

            % Regularization by shrinkage
            Sw = Sw.*(1-reg) + eye(nFeat).*reg;

            [W, LAMBDA] = eig(Sb,Sw);
            lambda = diag(LAMBDA);
            [lambda, SortOrder] = sort(lambda,'descend');

            % Handle ill conditioned data (nans)
            nNans = sum(isnan(lambda));
            SortOrder = circshift( SortOrder, -nNans);
            W = W(:,SortOrder);
            Y = X*W;



        end
        
        function clrs = GetColors( nClrs )
            try 
                clrs = cmap( nClrs );
            catch
                clrs = lines( nClrs );
            end
        end
        
        
    end
end
