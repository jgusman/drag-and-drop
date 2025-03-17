%% Feature selection analysis (Fig 5)
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
ddFigHelper.SetSaveDir(fullfile(saveFiguresFolder,'GH_FeatureSelection_out'))

%% Calculate Gesture and Attempt Selective Features and Top 400 Features Used in Latch Decoder (via MRMR)

kwpAllDurs = {};
selFeat = {};
topFeat = {};
numTopFeats = 400;
sigAlpha = 0.001;
selectivityTypes = {'gest','attempt'};

for ii = 1:length(sesData_GH_T11)
    [~, kwp_allDur{ii}] = CalcGestAndAttemptSelectiveFeatures(sesData_GH_T11(ii),selectivityTypes);
    [selFeat_T11{ii},topFeat_T11{ii}] = GetSelectiveAndTopFeats(sesData_GH_T11(ii),kwp_allDur{ii},sigAlpha,numTopFeats);
end

for ii = 1:length(sesData_GH_T5) %(should just be one session with T5)
    [~, kwp_allDur{ii}] = CalcGestAndAttemptSelectiveFeatures(sesData_GH_T5(ii),selectivityTypes);
    [selFeat_T5{ii},topFeat_T5{ii}] = GetSelectiveAndTopFeats(sesData_GH_T5(ii),kwp_allDur{ii},sigAlpha,numTopFeats);
end

%% Plot Figure 5(A):
ddFigHelper.LogPrints('Numbers of Selective Features by Feature Type')

disp('Significant Features - T11')
PlotFeatureProportionBarPlot(sesData_GH_T11,selFeat_T11);
ddFigHelper.SaveFigure('Figure 5A (T11) - Selective feats by feature type')

disp('Significant Features - T5')
PlotFeatureProportionBarPlot(sesData_GH_T5,selFeat_T5);
ddFigHelper.SaveFigure('Figure 5A (T5) - Selective feats by feature type')

ddFigHelper.LogPrints('off')


%% Plot Figure 5(B):
ddFigHelper.LogPrints('Numbers of Features Used in Decoding by Feature Type')

disp('Top Features - T11')
PlotFeatureProportionBarPlot(sesData_GH_T11,topFeat_T11,numTopFeats);
ddFigHelper.SaveFigure('Figure 5B (T11) - Feats used by feature type')

disp('Top Features - T5')
PlotFeatureProportionBarPlot(sesData_GH_T5,topFeat_T5,numTopFeats);
ddFigHelper.SaveFigure('Figure 5B (T5) - Feats used by feature type')

ddFigHelper.LogPrints('off')




%% %%%%% INTERNAL HELPER FUNCTIONS %%%%% %%


function PlotFeatureProportionBarPlot(sesData,selFeat,nTopFeats)

if nargin <3
    nTopFeats = [];
end

%% Plot Proportions of feats selective for gesture vs attempt vs both for each Feature Type
participant = upper(sesData(1).participant);

featPltOrder = [1 2	7 6	5 4	3];
nFeat = size(selFeat{1},1);
nFeatTypes = nFeat/192;

barColors = ddHelper.decoder.colors;
barColors(2,:) = barColors(3,:);

if isempty(nTopFeats)
    saveNameMod = 'selective';
else
    saveNameMod = ['Top ', num2str(nTopFeats)];
end

    n_gestSel_only = zeros(nFeatTypes,length(sesData));
    n_attemptSel_only = zeros(nFeatTypes,length(sesData));
    n_bothSel = zeros(nFeatTypes,length(sesData));

    for sesInd = 1:length(sesData)
        for fT = 1:nFeatTypes
            featInds = ((fT-1)*192+1) : (fT*192);
            
            selFeat_gest = selFeat{sesInd}(featInds,1);
            selFeat_attempt = selFeat{sesInd}(featInds,2);

            n_gestSel_only(fT,sesInd) = sum(and(selFeat_gest, ~selFeat_attempt));
            n_attemptSel_only(fT,sesInd) = sum(and(~selFeat_gest, selFeat_attempt));
            n_bothSel(fT,sesInd) = sum(and(selFeat_gest, selFeat_attempt));
        end

        % print results
        fprintf('\n%s %s Features (Session %d)\n\n', participant, upper(saveNameMod), sesInd)
        selectivityTable = table(n_gestSel_only(featPltOrder,sesInd),n_bothSel(featPltOrder,sesInd),n_attemptSel_only(featPltOrder,sesInd),...
            'VariableNames', {'Gesture Only','Both','Attempt Only'},...
            'RowNames',ddHelper.featTypeNames(featPltOrder));
        DisplayUnformatted(selectivityTable)
    end

    n_gestSel_only = mean(n_gestSel_only,2); %avg over sessions if needed (i.e for T11)
    n_attemptSel_only = mean(n_attemptSel_only,2);
    n_bothSel = mean(n_bothSel,2);
    
    % print results
    fprintf('\n%s %s Features (BOTH Sessions)\n\n', participant,upper(saveNameMod))
    selectivityTable = table(n_gestSel_only(featPltOrder,1),n_bothSel(featPltOrder,1),n_attemptSel_only(featPltOrder,1),...
        'VariableNames', {'Gesture Only','Both','Attempt Only'},...
        'RowNames',ddHelper.featTypeNames(featPltOrder));
    DisplayUnformatted(selectivityTable)

    barData = [n_gestSel_only(:),n_bothSel(:),n_attemptSel_only(:)];

    % Setup figure
    fwh = [3 1.5];
    ddFigHelper.CreateFigure([],fwh)
    clf
    axArgs = ddFigHelper.GetDefaultAxArgs();
    axArgs.marg_w = [0.13 0.05]; % Update axes width boundaries
    axArgs.marg_h = [0.2 0.08]; % Update axes height boundaries
    axs = ttight_subplot(1,1,axArgs);
    ax = axs;

    barData = barData(featPltOrder,:);
    bchart = bar(barData,'stacked','BarWidth',0.7);

    for bb = 1:length(bchart)
        bchart(bb).FaceColor = barColors(bb,:);
        bchart(bb).EdgeColor = barColors(bb,:);
        bchart(bb).LineWidth = 0.5;
    end
    bchart(2).EdgeColor = barColors(1,:);
    hF = hatchfill2(bchart(2),'single','HatchAngle',45,'HatchLineWidth',1,'HatchColor',barColors(1,:)); %draw hatch lines over "both" condition

    xticklabels(ddHelper.featTypeNames(featPltOrder))
    xticks(1:nFeatTypes)
    ax.XAxis.FontSize = 8;
    ylabel('# of feats')
    xlim([0.4 nFeatTypes+0.6])
    ylim([0 192])
    yticks([0 96 192])

%     legend  %legend doesnt work with the hatch. make in Illustrator
%     uistack(bchart(1),"bottom")
%     set(axs,'Children',[bchart(3),hF,bchart(2),bchart(1)])


end


function DisplayUnformatted(stuff,title)   

   if nargin > 1
       fprintf('%s\n\n',title); 
   end

   stuffToDisp = formattedDisplayText(stuff,"SuppressMarkup",true);

   disp(stuffToDisp);

end