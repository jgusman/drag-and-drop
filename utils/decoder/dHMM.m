classdef dHMM < matlab.mixin.Copyable
    % This is an HMM decoder
    % Accessible functions:
    %
    % Train(obj, trainInds)
    % [predState, liklihoods] = Test(obj, testData, testLabels)
    % likelihoods = Predict(obj, testData)
    % likelihoods = Viterbi(obj,data)
    % xVal(obj)
    % Summarize(obj)
    properties
        data  % neural data, [nSteps x nFeat]
        labels % discrete labels for each class,  [nSteps x 1]
        
        eventInds
        performance % struct containing performance results
        coefficients % struct with parameters, coefficients for the decoder
        params % Defines how to calc the multistate decoder        
        info
    end
    
    
    
    methods
        function obj = dHMM( params )
        
            % Default parameters            
            defaultParams.initProbs                    = []; % Set init probabilities
            defaultParams.pureHMM                      = 0; % If true, compute hmm through em
            defaultParams.eventInds                    = [];
            defaultParams.enforceTransitionMatrix      = []; % If empty, uses empirical trans matrix. If scalar, sets trans matrix diag to scalar value and other 1-scalar. If matrix, sets to matrix.
            defaultParams.smoothBeforeNormalizing      = false;
            defaultParams.smoothLikFun                 = @movmean; %'ExpSmooth'; % @movmean; % movmean, GaussSmooth2. For smoothing the likelihoods
            defaultParams.smoothLikAmount              = [0 0]; % Additional param to smoothLikFun [nBack nForward]
            defaultParams.contrainStateCov             = ''; % Can be empty, 'diag', 'same'
            defaultParams.minProb                      = [];
			defaultParams.numStates                    = []; % If set, we force transition matrix, mu, emission to be this size
            if nargin < 1
                params = [];
            end
            % Init event inds
            obj.eventInds.train = [];
            obj.eventInds.test = [];
            
            obj.params = MergeBintoA(defaultParams,params);
            

%             if isempty(which('mixgauss_init'))
%                 error('HMM Toolbox not on path! This decoder requires Kevin Murphy''s HMM Toolbox.')
%             end
        end
        
        
        
        function Train(obj, trainInds)
        % Calibrate the Multistate decoder
            if nargin < 2 
                trainInds = [];
            end
            if ~isempty(trainInds)
                obj.eventInds.train = trainInds;
            else
                obj.eventInds.train = obj.params.eventInds;
            end
            
            % If no train inds, set to length of data
            if isempty(obj.eventInds.train)
                obj.eventInds.train = 1:size(obj.data,1);
            end
            

            % Compute HMM params
            [eMu, eCov, transitionMatrix] = CalculateStateMeansVarsAndTransitionMatrix(obj);
            
            
            obj.info.train.empiricalTransitionMatrix = transitionMatrix; 
            obj.coefficients.transitionMatrix = transitionMatrix;
            obj.coefficients.transitionMatrix = ReturnTransitionMatrix(obj);
            obj.coefficients.emissionMu       = eMu;
            obj.coefficients.emissionCov      = eCov;            
            
%             stateProb = obj.ComputeStateLWithToolbox(projData, obj.coefficients.emissionMu, obj.coefficients.emissionCov, obj);
        end
        
        
        function [predState, liklihoods, obsLiks] = Test(obj, testData, testLabels)
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
            testInds = obj.eventInds.test;
            if isempty(testInds), testInds = 1:size(obj.data,1); end
            testData = obj.data(testInds,:);
            testLabels = obj.labels(testInds);
        end
        

        % Make sure test labels match train label order
        if ~iscategorical(testLabels)
            testLabels = categorical(testLabels);
        end
        testLabels = MakeCatCommon(obj.labels, testLabels);
%         testLabels = removecats(testLabels);
        
        [hmmLiks,obsLiks] = ComputePred( obj, testData );
        
        % Metrics computed and saved to obj.performance
        [perf, info] = obj.ComputeMetrics( hmmLiks, testLabels );
        
        % Save off for summary
        obj.performance              = perf;
        obj.info.test.labels         = info.labels;
        obj.info.test.uStates        = info.uStates;
        obj.info.test.estLiks        = info.estLiks;
        obj.info.test.obsLiks        = obsLiks;
        obj.info.test.predState      = info.predState; 
        obj.info.test.incorrectGuess = info.incorrectGuess;
        
        if nargout > 0
            predState  = info.predState;
            liklihoods = info.estLiks;
        end
        end
        
        function [hmmLiks,obsLiks] = Predict(obj, testData)
        
            if nargin < 2
                testData = obj.data(1:size(obj.data,1),:);
            end

            [hmmLiks,obsLiks] = ComputePred( obj, testData );        
        
        end
        
        function likelihoods = Viterbi(obj,data)
            likelihoods = ComputeViterbi(obj, data);
        end
        
        function xVal(obj)
        
            if ~isempty(obj.params.kFoldIter)
                kFoldIter = obj.params.kFoldIter;
            else
                kFoldIter = obj.params.kFold;
            end
            
            metric  = nan(kFoldIter,1);
            
            
            % Perform cross validation
            iter = 1;
            done = 0;
            while ~done
                
                cvInds = GetFoldInds(obj);
                
                for kk = 1:obj.params.kFold
                    
                    testInds  = cvInds == kk;
                    trainInds = find(~testInds);

                    obj.Train(trainInds);

                    obj.Test( obj.data(testInds,:), obj.labels(testInds,:) );

                    metric(iter) = obj.performance.auc(2);
                    xPerformance(iter) = obj.performance;
                    
                    iter = iter + 1;
                    done = iter > kFoldIter;
                    if done
                        break
                    end

                end
            end
            
            % Average the performance
            % For each performance metric, concatenate the metric then
            % average.
            fn = fieldnames(obj.performance);
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
        
        function cvInds = GetFoldInds(obj)
            nPoints = size(obj.data,1);
            
            
            if ~isempty(obj.params.trialBoundaries)
                if isvector(obj.params.trialBoundaries)
                    % assume trialBoundaries are the start inds
                    % final index is the equivalent to the stop trial index
                    % + 1
                    tb = obj.params.trialBoundaries;
                    for ii = 1:length(obj.params.trialBoundaries)-1
                        trialStartStops(ii,:) = [tb(ii), tb(ii+1)-1];
                    end
                    
                else
                    trialStartStops = obj.params.trialBoundaries;
                end
                
                
                nTrials = size(trialStartStops,1);
                kFold = obj.params.kFold;
                % Divide the trials into the folds
%                 if obj.params.crossValEqualStatePerFold
%                 else
%                 end
                %%
                
                lblPerTrial = obj.labels(trialStartStops(:,1),:);
                uLbls = unique(lblPerTrial);
                cvTrls = nan(nTrials,1);
                for ii = 1:length(uLbls)
                    lblTrlsInds = find(lblPerTrial == uLbls(ii));
                    nTrialInds = length(lblTrlsInds);
                    nPerFold = floor( nTrialInds/kFold );
                    if nPerFold == 0
                        % Not a full set, set the first N folds
                        cvTrls(lblTrlsInds) = 1:nTrialInds;
                    else
                        for jj = 1:kFold
                            foldInds = lblTrlsInds( randperm(length(lblTrlsInds), nPerFold) ); 
                            lblTrlsInds = setdiff(lblTrlsInds,foldInds);
                            cvTrls(foldInds) = jj;
                        end
                        % Any remaining inds go into the last fold
                        cvTrls(lblTrlsInds) = jj;
                    end
                    
                end
                %%
%                 cvTrls = crossvalind('Kfold',nTrials,kFold);
                
                %% Set each trial to it's associated fold number
                cvInds = nan(nPoints,1);
                for ii = 1:length(cvTrls)
                    trialInds = trialStartStops(ii,1):trialStartStops(ii,2);
                    cvInds(trialInds) = cvTrls(ii);
                end
            
            else
                cvInds = crossvalind('Kfold',nPoints,obj.params.kFold);
            end
        end
        
        function Summarize(obj)
        % Plot/summarize/save results
        
        
        end
        
        function sobj = saveobj(obj)
        % Save property fields, but do not save the object.
        % TO DO: have a function that can load this struct to recreate the
        % Multistate object
             fn = fieldnames(obj);
             for ii = 1:length(fn)
                 sobj.(fn) = obj.(fn);
             end

             if ~obj.params.onSave.keepDataAndLabels
                 sobj.data = [];
                 sobj.labels = [];
             end    
        end
        
        
    end
    
    methods (Access = protected)
        
        function transitionMatrix = ReturnTransitionMatrix(obj)
            % Use an enforced transition matrix
            %    if params.enforceTransitionMatrix is not empty
            % Use the empirical computed transition matrix
            %   if it is not empty and enforced transition matrix is not empty        
            % Use the current contents of coefficients if both are empty
            try
                empiricalTransitionMatrix = obj.info.train.empiricalTransitionMatrix;
            catch
                empiricalTransitionMatrix = [];
            end

            % Default if everything else is empty
            transitionMatrix = obj.coefficients.transitionMatrix;

            % Look to see if obj.params.enforceTransitionMatrix is empty
            % If not, use the empiricalTransitionMatrix
            if isempty(obj.params.enforceTransitionMatrix)
                if ~isempty(empiricalTransitionMatrix)
                    transitionMatrix = empiricalTransitionMatrix;
                end
            else
                % Use the param enforceTransitionMatrix to set the transition
                % matrix
                nStates = size(transitionMatrix,1);
                transitionMatrix = zeros(nStates);

                % Scalar value, transition matrix = I * scalar
                if isscalar( obj.params.enforceTransitionMatrix )
                    if obj.params.enforceTransitionMatrix > 0
                        
%                         diagInds = logical(eye(nStates));
%                         transitionMatrix(diagInds) = obj.params.enforceTransitionMatrix;
%                         transitionMatrix(~diagInds) = (1-obj.params.enforceTransitionMatrix)/(nStates-1);
                        
                        % Get number of trained states
                        if ~isempty(obj.labels)
                        
                            obj.labels = categorical(obj.labels);
                            uStates  = categories( obj.labels );
                            nTrainedStates = length(uStates);  % construct enforced tran matrix only using relevant states (because often forced to use 11 states for model)
                        elseif ~isempty(empiricalTransitionMatrix)
                            nTrainedStates = sum(empiricalTransitionMatrix(1,:) > 0);
                        else
                            nTrainedStates = nStates; %just use nStates if can't find specific trained states
                            warning("Couldn't find nTrainedStates, so using nStates. Input complete enforced transition matrix to ensure expected trans matrix used.")
                        end
                        
                        % Create transmat based on nTrainedStates
                        transMat_temp = zeros(nTrainedStates);
                        diagInds = logical(eye(nTrainedStates));
                        transMat_temp(diagInds) = obj.params.enforceTransitionMatrix;
                        if nTrainedStates > 1
                            transMat_temp(~diagInds) = (1-obj.params.enforceTransitionMatrix)/(nTrainedStates-1);
                        end
                        % Update Values to full transition matrix
                        transitionMatrix(1:nTrainedStates,1:nTrainedStates) = transMat_temp;

                    else
                        % Zero the transition matrix
                    end
                else
                    % Set transitionMatrix to enforceTransitionMatrix
                    if any(~ismember( nStates, size(obj.params.enforceTransitionMatrix) ))
                        error(sprintf('Enforced transition matrix is not the same size as nStates.\nnStates: %d, enforcedTransitionMatrix: %s\n', ...
                            nStates, num2str(size(obj.params.enforceTransitionMatrix))))
                    end
                    transitionMatrix = obj.params.enforceTransitionMatrix;
                end

            end
        end
        
        function mostLikelihood = ComputeViterbi(obj, data )
            
            [~, ~, trans] = GetMeansVarsAndTransitionMatrix(obj);
            
            prior = obj.params.initProbs;
            if isempty(prior)
                nStates = size(trans,1);
                prior= zeros(1,nStates);
                prior(1) = 1;
            end
            
            likelihoods = ComputePred( obj, data );
            mostLikelihood = viterbi_path(prior, trans, likelihoods(:,1:nStates)')';
        end
        
        function [hmmLiks,obsLiks] = ComputePred( obj, testData )
            initProbs = obj.params.initProbs;
            
            nPoints = size(testData,1);
            nStates = size(obj.coefficients.emissionMu,1);
            
            useToolbox = 0;
            if useToolbox
                stateProb = obj.ComputeStateLWithToolbox(testData, obj.coefficients.emissionMu, obj.coefficients.emissionCov, obj);
            else
                
                if isempty(initProbs)
                    initProbs= zeros(1,nStates);
                    initProbs(1) = 1;
                end

                hmmLiks = zeros(nPoints,nStates);
                %% Toolbox prediction method (much faster)
                [mu1, Sigma1, trans, ignoreStateDims, ignoreOutDims] = GetMeansVarsAndTransitionMatrix(obj);
                
                testData(:,ignoreStateDims) = [];
                initProbs(ignoreOutDims) = [];
                obsLiks = mixgauss_prob(testData', mu1, Sigma1);

                smoothFun = GetSmoothFun(obj);
                if obj.params.smoothBeforeNormalizing
                    obsLiks = smoothFun(obsLiks',obj.params.smoothLikAmount)';
                end
                
                if all(trans(:)==0)
                    tmpProbs = normalise(obsLiks,1);
                else
                    tmpProbs  = fwdback(initProbs, trans, obsLiks,...
                                    'fwd_only', 1, 'scaled', 1, 'handle_zero_lik', 1);
                end
                
                hmmLiks(:,~ignoreOutDims) = tmpProbs';
                
                %% Straightforward, analogous to model
%                 tt = tic;
%                 probs = obj.RunForwardModel( d, ...
%                                              initProbs, ...
%                                              ReturnTransitionMatrix(obj), ...
%                                              obj.coefficients.emissionMu, ...
%                                              obj.coefficients.emissionCov);
%                 toc(tt)
%                 Compare
%                 figure
% %                 plot( probs(:,st) )
% %                 plot( stateProb(:,st) )
%                 plot( sum(probs- stateProb,2) )
                

                %% post estimate smoothing
                if ~obj.params.smoothBeforeNormalizing
                    hmmLiks = smoothFun( hmmLiks, obj.params.smoothLikAmount );
                end
                
                % Keep consistent with hmmLiks - nSteps x nStates
                obsLiks = obsLiks';
                

            end
        end
        
        function smoothFun = GetSmoothFun(obj)
            if ischar(obj.params.smoothLikFun)
                if ismember( obj.params.smoothLikFun, methods(obj) )
                    % Assume we are calling the local smooth function
                    smoothFun = @(varargin) obj.(obj.params.smoothLikFun)(varargin{:});
                else
                    % Assume it is a regular function
                    smoothFun = str2func( obj.params.smoothLikFun );
                end
            else
                smoothFun = obj.params.smoothLikFun;
            end
        end
        
        function [mu, sigma, trans, ignoreStateDims, ignoreOutDims] = GetMeansVarsAndTransitionMatrix(obj)
            % Returns transposed from our normal set up such that it is
            % ready to be used by the HMM toolbox.
            
            mu = obj.coefficients.emissionMu';
            ignoreStateDims = all(mu==0,2) | any(isnan(mu),2);
            ignoreOutDims = all(mu==0,1) | any(isnan(mu),1);
            mu(:,ignoreOutDims) = [];
            mu(ignoreStateDims,:) = [];

            sigma = permute(obj.coefficients.emissionCov, [2 3 1]);
            sigma(:,:,ignoreOutDims) = [];
            sigma(ignoreStateDims,:,:) = [];
            sigma(:,ignoreStateDims,:) = [];

            

            trans = ReturnTransitionMatrix(obj);
            trans = trans(~ignoreOutDims,~ignoreOutDims);
            
        end
        
        function [eMu, eCov, transitionMatrix] = CalculateStateMeansVarsAndTransitionMatrix(obj)
            obj.labels = categorical(obj.labels);
            uStates  = categories( obj.labels );
            grpInd   = double(obj.labels);
            projData = obj.data;
            if ~isempty(obj.params.numStates)
                nStates = obj.params.numStates;
            else
                nStates = length(uStates);
            end

            
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
            nComp = size(obj.data,2);
            eMu   = nan(nStates, nComp);
            eCov  = nan(nStates, nComp, nComp);
            
            for ii = 1:nStates
                eMu(ii,:)    = nanmean(projData(grpInd == ii,:));
                eCov(ii,:,:) = nancov(projData(grpInd == ii,:));
            end
            
            % Covariance constraints
            if ~isempty(obj.params.contrainStateCov)
                switch obj.params.contrainStateCov
                    case 'diag'
                        for ii = 1:nStates
                            eCov(ii,:,:) = diag(diag(squeeze(eCov(ii,:,:))));
                        end
                    case 'same'
                        sameCov = mean(eCov,1);
                        for ii = 1:nStates
                            eCov(ii,:,:) = sameCov;
                        end
                    otherwise
                        error('Unrecognized Cov constraint, %s', obj.params.contrainStateCov) 
                end
            end


        end
    end
    
    methods (Static)
        function [perf, info] = ComputeMetrics(estLiks, testLabels, ignoreFirstState)
        % Computes the following metrics (obj.performance)
        % auc
        % confustion matrix
        % precision
        % recall
        % F-measure
       
       
            if nargin < 3
                ignoreFirstState = false;
            end
       
       
            testLabels = categorical( testLabels );
            grpInd = double(testLabels(:));

            nPredStates = size(estLiks,2);
            nPoints = length(grpInd);
            nStates = length(unique(grpInd(~isnan(grpInd))));
            nStates = max(nStates, nPredStates);
        
            
            rocThresh = [];
            auc = [];
            for ii = 1:nPredStates
                nClassObs = sum(grpInd==ii);
                if nClassObs~=0 && nClassObs ~= length(grpInd)
                    [X,Y,T,AUC, optro] = perfcurve(grpInd==ii,estLiks(:,ii), 1);
                    rocThresh(ii) = T( X==optro(1) & Y==optro(2) );
                    auc(ii) = AUC;
                else
                    rocThresh(ii) = nan;
                    auc(ii) = nan;
                end
            end

            %%
            
            [maxProb,predState] = max(estLiks,[],2);
            
%             thresh = 0.95;
%             predState(maxProb<thresh) = 1; % Assume first state is null

            incorrectGuess = double(grpInd(:)') ~= double(predState(:)');
            
            
            % Save off for summary
            info.labels    = testLabels;
            info.uStates   = categories(testLabels);
            info.estLiks   = estLiks;
            info.predState = categorical( predState, 1:length(info.uStates), info.uStates );
            info.incorrectGuess = incorrectGuess;
            
            
            confusionMatrix = confusionmat(info.labels,info.predState);
            %% Compute metrics           
            try
            fn = zeros(nStates,1);
            fp = zeros(nStates,1);
            tn = zeros(nStates,1);
            tp = diag(confusionMatrix);
            for ii = 1:nStates
                negInds                    = true(nStates,1);
                negInds(ii)                = false;
                tnLogical                  = false(nStates);
                tnLogical(negInds,negInds) = true;
                
                fn(ii) = sum( confusionMatrix(ii,negInds) );
                fp(ii) = sum( confusionMatrix(negInds,ii) );
                tn(ii) = sum( confusionMatrix(tnLogical), 'all' );
            end
            catch err
                UnrollError(err)
                keyboard
            end
          

%%
            recall = tp./(tp+fn); 
            recall(isnan(recall)) = 0;
            precision = tp./(tp+fp);
            precision(isnan(precision)) = 0;
            accuracy = (tp + tn)./(tp+tn+fp+fn);
            fmeas = 2*(precision.*recall)./(precision+recall);
            fmeas(isnan(fmeas)) = 0;
            
            %% MCC - Matthews Correlation Coefficient
            % Comparing two K-category assignments by a K-category correlation coefficient
%             mcc = ComputeMCC(conf);
            % Numerator
            tmpSum = 0;
            for kk = 1:nStates
                for ll = 1:nStates
                    tmp = confusionMatrix(kk,:)*confusionMatrix(:,ll);
                    tmpSum = tmpSum + tmp;
                end
            end
            numerator = nPoints*trace(confusionMatrix) - tmpSum;
            
            % Denominator
            tmpSum1 = 0;
            tmpSum2 = 0;
            for kk = 1:nStates
                for ll = 1:nStates
                    tmpSum1 = tmpSum1 + confusionMatrix(kk,:)  * confusionMatrix(ll,:)';
                    tmpSum2 = tmpSum2 + confusionMatrix(:,kk)' * confusionMatrix(:,ll);
                end
            end
            den1 = sqrt( nPoints^2 - tmpSum1 );
            den2 = sqrt( nPoints^2 - tmpSum2 );
            
            mcc = numerator / ( den1 * den2 );

            
            %% Set performance values
            perf.rocThresh = rocThresh;
            
            perf.auc = auc;
            perf.conf = confusionMatrix;    

            perf.tn = tn;
            perf.tp = tp;
            perf.fn = fn;
            perf.fp = fp;
            
            perf.accuracy = accuracy;
            perf.recall = recall;
            perf.precision = precision;
            perf.fmeas = fmeas;
            perf.mcc = mcc;
        end

        function probs = RunForwardModel( observation,initialStateProb, transitionMatrix, emissionMu, emissionCov)
        % Forward Algorithm
        %
        % USAGE: out = hmmForwardAlgorithm(observation,initialStateProb,transitionMatrix,emission)
        %
        % Inputs
        %       observation      - [nPoints x nStates]
        %       initialStateProb - [1 x nStates] initial Probabilities for each state
        %       transitionMatrix - [nStates x nStates]
        %       emissionMu       - [nState x nDim]
        %       emissionCov      - [nStates x nDim x nDim]
        %
        % Output
        %       probs            - [nPoints x nStates] pred probs
        %
        %%

            nPoints  = size(observation,1);
            nStates = numel(initialStateProb);
            dims    = find( sum(abs(emissionMu)) ~= 0 );
            probs(1,1:nStates) = initialStateProb;%[Pclick Pmove]
            obserProbGstate    = zeros(1,nStates);

            for t = 2:nPoints
                for iState = 1:nStates
                    c = squeeze( emissionCov(iState,dims,dims) );
                    if all(sum(c)==0)
                        continue;
                    end
                    m = squeeze( emissionMu(iState,dims) );
                    obserProbGstate(iState) = mvnpdf(observation(t,dims),m,c);
                end

                probs(t,:) = (probs(t-1,:)*transitionMatrix).*obserProbGstate;
                probs(t,:) = probs(t,:)./sum(probs(t,:));

                if isnan(probs(t,:))
                    keyboard;
                end
            end
        end
        function stateProb = ComputeStateLWithToolbox(projFeatures, emissionMu, emissionCov, obj)
        
          O = size(projFeatures,2);          %Number of coefficients in a vector 
          T = size(projFeatures,1);         %Number of vectors in a sequence 
          nex = 1;        %Number of sequences 
          M = 1;          %Number of mixtures 
          Q = size(emissionMu,1);          %Number of states 

        % [guessTR,guessE,logliks] = hmmtrain(grpInd,transitionMatrix,em);
        % projFeatures
        data = projFeatures';
        % zz data
        % data = randn(O,T,nex);
        %%
        [mu0, Sigma0] = mixgauss_init(Q*M, data, 'full');
        mu0 = reshape(mu0, [O Q M]);
        Sigma0 = reshape(Sigma0, [O O Q M]);
        mixmat0 = mk_stochastic(rand(Q,M));
        prior0 = normalise(rand(Q,1));
        transmat0 = mk_stochastic(rand(Q,Q));
        %%
%         mu0(:,:,1) = emissionMu';
%         Sigma0(:,:,:,1) = permute(emissionCov, [2 3 1]);
        [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
            mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', 50);
        % 

%%
        [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
            mhmm_em(data, prior0, transmat0,emissionMu', permute(emissionCov, [2 3 1]), mixmat0, 'max_iter', 50, ...
            'adj_mu', 0, 'adj_Sigma', 0);



        %% Pred only
        d = data; %(:,1:10e3);
        mu1 = emissionMu';
        Sigma1 = permute(emissionCov, [2 3 1]);
        mixmat1 = mixmat0;
        
        prior1 = prior0;
        tFull = tic;
        obslik = mixgauss_prob(d, mu1, Sigma1, mixmat1);
        mixTime = toc(tFull)
        [alpha, beta, gamma, ll] = fwdback(prior1, transmat1, obslik, 'fwd_only', 1);
        fullTime = toc(tFull)
        
        %%
        d = data;
%         d = obj.data(obj.eventInds.test,:) * obj.coefficients.projection;
%         d = d';
        obslik = mixgauss_prob(d, mu1, Sigma1, mixmat1);
        [alpha, beta, gamma, ll] = fwdback(prior1, transmat1, obslik, 'fwd_only', 1);


        stateProb = alpha';
        %%
        [~, predState] = max(stateProb');
%         labels = obj.labels(obj.eventInds.train);
        labels = obj.labels(obj.eventInds.test);
        sum(labels(:)==predState(:))./length(labels)
        %%
        figure
        ConfMat = confusionchart(labels(:)',predState(:)', 'RowSummary','total-normalized');
        %%
        figure
        plot(labels)
        plot( stateProb )
        end
    end
end


