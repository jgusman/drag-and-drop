function [repeatWindowSize, windowStartStops] = FindConsecutiveOnes( vectVar, includeSingleInstances )
% function [repeatWindowSize, windowStarts] = FindConsecutiveOnes( logicalVect, includeSingleInstances )
% 
% This function find the windows of consecutive non-zero values, and
% outputs the length of the window as well as where the window starts
% 
% Inputs:
%   vectVar 
%       A 1D vector of values that you want to count non-zero windows
% 
%   includeSingleInstances (default 0)
%       Flag, if true, single instances of non-zero will be included in the
%       output.
% 
% Output:
% 
%   repeatWindowSize
%       An n x 1 vector where each row corresponds to the window length of
%       the windowStartStops.
% 
%   windowStartStops
%       An n x 2 matrix of start, stop indices of non-zero windows inside
%       vectVar



if nargin < 2
    includeSingleInstances = 0;
end
plotOn = 0;


oneInds                 = find(vectVar(:)');                                % Where are the non-zero values
repeatSeparation        = diff(oneInds);                                    % How far away are the non-zero values from each other
consecutiveInds         = [false repeatSeparation == 1 false];              % Diff == 1, i.e. repeated non-zero values
consecutiveWinInds      = find( diff(consecutiveInds) );
consecutiveWindowPairs  = reshape(consecutiveWinInds,2,[]);                 % Start stop of the repeats


if includeSingleInstances && ~isempty(oneInds)

    %%
    % Find single instances.
    singleInstanceInds = find(~consecutiveInds);
    singleInstanceInds = singleInstanceInds(diff(singleInstanceInds)==1);
    singleWinPairs     = repmat(singleInstanceInds,2,1);
    
    % Combine with consecutive inds
    consecutiveWindowPairs = cat(2, singleWinPairs, consecutiveWindowPairs);
    
    % Sort
    [~, sortI]             = sort( consecutiveWindowPairs(1,:) );
    consecutiveWindowPairs = consecutiveWindowPairs(:,sortI);
end


% Reshape for output
repeatWindowSize        = diff(consecutiveWindowPairs) + 1;                 % Get the size of each window (diff between the starts and stops)
windowStartStops        = oneInds(consecutiveWindowPairs)';                 % Output index start stop pairs


% In the scenario where there is just one set of consecutive ones,
% transpose the window.
if size(windowStartStops,2) == 1
    windowStartStops = windowStartStops';
end


%% Plot
if plotOn
    plotWin = 10;
    PlotWindow(vectVar, windowStartStops, repeatWindowSize, plotWin)
end
end

function PlotWindow(upsampleStepsLogical, windowStartStops, repeatWindowSize, plotWin)


figure;
a = find(repeatWindowSize > plotWin, 1);
if ~isempty(a)
    winSize = repeatWindowSize(a);
    logicalIndex = windowStartStops(a,1);
    pltInds = logicalIndex-1:logicalIndex+winSize;
    figure;
    plot(upsampleStepsLogical(pltInds))
    title(sprintf('repeat window length: %d', winSize ));
    sumWindow = sum(upsampleStepsLogical(pltInds));
    fprintf('Window size found %d, indexed %d. Is equal: %d\n', winSize, sumWindow, isequal( winSize, sumWindow ));
else
    fprintf('Could not find a consecutive window > %d \n', plotWin);
end
end