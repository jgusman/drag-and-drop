function labels = GetLabelsForStartStops(varargin)
%% Repeat trial labels given start stops
% labels = GetLabelsForStartStops(taskInfo) 
% OR
% labels = GetLabelsForStartStops(labelsPerTrial,startStops)
% 
% taskInfo.trialStartStops
% taskInfo.labelPerTrial


if length(varargin) == 1 && isstruct(varargin{1})
    taskInfo = varargin{1};
    if isfield(taskInfo, 'labelPerTrial')
        labelsPerTrial = taskInfo.labelPerTrial;
    elseif isfield(taskInfo, 'labels')
        labelsPerTrial = taskInfo.labels;
    else
        error('Label field was not found!')
    end
    
    if isfield(taskInfo, 'trialStartStops')
        startStops = taskInfo.trialStartStops;
    elseif isfield(taskInfo, 'startStops')
        startStops = taskInfo.startStops;
    else
        error('StartStop field was not found!')
    end
    
elseif length(varargin) == 2
    labelsPerTrial = varargin{1};
    startStops = varargin{2};
else
    error('Unknown input type')
end


% Set labels to have the same category set and order as labelPerTrial
catStr   = categories(labelsPerTrial);
labels   = categorical([],1:length(catStr), catStr);
ststpLen = diff(startStops,[],2)+1;

% Concat repeated label for each trial's start stop inds
for ii = 1:length(labelsPerTrial)
    repLabels = repmat( labelsPerTrial(ii), ststpLen(ii), 1 );
    labels    = cat(1, labels, repLabels);
end

end