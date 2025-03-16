classdef ddHelper
% ddHelper
% 
% Static helper class for drag and drop analyses.

    properties (Constant)
        gestureHeroDates = {'2020.07.21', '2020.07.24'}; % Booted '2020.04.10'
        gestureHeroDatesT5 = {'2023.04.04'}; 
        dragAndDropDates = {'2021.03.04', '2021.03.09'}; % Booted '2021.02.16',

        
        excludeNs5OutlierTrialsForPerformance = true;
        
        gestureColors = cat(1,ones(1,3).*0.3,...
                              lines(8));
                          
        decoder = struct('colors',[107, 133, 183
                                   132, 169, 99
                                   198, 149, 52]./255, ...
                         'names', {{'Gesture', 'Latch', 'Attempt'}});

        % Used for gesture hero durations
        % 1, 2, 4 second holds
%         durationColors = [0.2200    0.6364    0.9892
%                           0.1050    0.4694    0.7781
%                           0.0600    0.2682    0.4446];

         durationColors =      [88 154 211
                                39 120 188
                                16 69 114]./255;

         attemptDurationColors =    [231 192 108
                                     227 156 39
                                     188 130 43]./255;
         trialTypeColors = [219 185 0           % move only
                            219 56 0            % click
                            219 133 0]./255;    % drag

         ddSesColors = [0.5 1 0.5
                         0.5 0.5 1];

         featTypeColors_attempt =   [233, 208, 34;
                                     233, 175, 30;
                                     232, 142, 26;
                                     232, 110, 22;
                                     231, 77, 17;
                                     231, 44, 13;
                                     230, 11, 9]./255;
         featTypeColors_gest =   [67, 193, 151;
                                     61, 164, 140;
                                     54, 136, 129;
                                     48, 107, 118;
                                     41, 78, 106;
                                     35, 50, 95;
                                     28, 21, 84]./255;
         featTypeColors =            [0, 0, 0;
                                     125, 125, 125;
                                     190, 184, 76;
                                     130, 189, 84;
                                     40, 133, 52;
                                     19, 89, 91;
                                     57, 149, 184]./255;


        
         featTypeNames = {'TX','SP','0-11Hz','12-19Hz','20-39Hz','40-128Hz','129-250Hz'};        

                   
    end
    
    methods (Static)
        function sesTaskInfos = GetTaskInfoForSessions( sessionInfos, featNames )
            % sesTaskInfos = GetTaskInfoForSessions( sessionInfos, featNames )
            % Gets session info (task, slc) excluding features.
            % sessionInfos is from GetDragAndDropSessionInfo
            
            if nargin < 2
                featNames = {};
            end
            
            %%
            
            for sesI = 1:length(sessionInfos)
            try
                sesInfo = sessionInfos(sesI);
                
                s = SummarizeSession( sesInfo.path );
                blks = sesInfo.blocks;
                if ~isempty(featNames)
                    [feat,slc] = ConcatNeuralData(s,blks,featNames);
                    [taskInfo,~,gameData] = sesInfo.taskFunction(s,blks,slc);
                else
                    [taskInfo,slc,gameData] = sesInfo.taskFunction(s,blks);
                end
                
                
                sesInfo.s = s;
                sesInfo.taskInfo = taskInfo;
                sesInfo.slc = slc;
                sesInfo.gameData = gameData;
                
                if ~isempty(featNames)
                    sesInfo.feat = feat;
                    sesInfo.featInfo = GetFeatInfo(slc,featNames);
                end
                
                sesTaskInfos(sesI) = sesInfo;
            catch err
                UnrollError(err)
                beep
                fprintf(2,'Pausing at error...\n');
                keyboard
            end
            end
        end
        
        
        function excludeTrials = ExcludePerformanceTrialsWithNS5Outliers(taskInfo)
            threshold = 5;
            if ddHelper.excludeNs5OutlierTrialsForPerformance
                % exclude epochs with > threshold
                excludeTrials = taskInfo.prctOutliersPerEpoch > threshold;
            else
                % No exclusion
                excludeTrials = false(size(taskInfo.prctOutliersPerEpoch));
            end
        end
    end
    
    
    
    %% Drag and drop info
    
    methods (Static)
        function PrintNumberOfTrialTypesPerBlock(taskInfo)
            %%
            uBlks = unique(taskInfo.blockNumber);
            uTrialTypes = categories(taskInfo.trialAttemptType);
            uTrialTypes(end) = [];
            nPad = 5;
            uTrialTypeStr = pad(uTrialTypes,nPad,'both');

            ii = 1;
            blockStr = sprintf('%02d) Block %02d ', ii,uBlks(ii));
            fprintf('%s %s %s %s\n', pad('Block #',length(blockStr),'both'), uTrialTypeStr{:})
            for ii = 1:length(uBlks)
                selI = ismember(taskInfo.blockNumber, uBlks(ii));
                selI = selI & contains(taskInfo.trialStage,'Prepare');
                nTrialsStr = {};
                for jj = 1:length(uTrialTypes)
                    selJ = selI & ismember(taskInfo.trialAttemptType,uTrialTypes(jj));
                    nTrialsStr{jj} = pad(num2str( sum(selJ) ),nPad,'both');
                end
                blockStr = sprintf('%02d) Block %02d ', ii,uBlks(ii));
                fprintf('%s %s %s %s\n', blockStr, nTrialsStr{:});
            end
        end
        
        function PlotGestureLegend(taskInfo)
            
        end
    end
end
