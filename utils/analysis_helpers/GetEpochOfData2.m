function varargout = GetEpochOfData2(data, events, window, mode, removeOutOfBounds)
% [epochData, epoch] = GetEpochOfData2(data, events, window, mode, removeOutOfBounds)
% 
% Given event inds and a window, indexes into data and outputs in standard
% shapes.
% 
% Data can be a cell array, e.g., {data, labels}. If that is the case, then
% there is an output for each element, and the last output is always epoch:
% [epochData, epochLabels, epoch] = GetEpochOfData({data, labels}, events, window, mode)
% 
% Inputs:
% 
%  data [nPoints x nDim] OR [nEvents x nDim] OR cell array of elements.
%    data to be indexed / reshaped.
% 
%  events
%    Vector of events to align the window vector, window. E.g. trial start
%    inds.
% 
%  window
%    Vector of offset inds from events used to index into data, e.g. 0:49
% 
%  mode (default: 'td')
%    Char or cell array (for each input) specifying the output shape.
% 
%    How the output is reshaped using the following letters:
%    t = Time, e = Event, d = Dim
% 
%    Examples
%    'td'  % Time-event x Dim
%    'ted' % Time x Event x Dim
%    'ed'  % Event x Dim
% 
%    For two element input, you can specify output shape as a cell array
%    output = {'ted', 'ed'};
%
%    If data is [nEvents x nDim] and output has a time component ('td', ted')
%    then, we will repeat data for each time instance.
% 
%  removeOutOfBounds (default true) 
%    If true, then events will be removed when events + window is out of
%    bounds to data (< 1 or > size(data,1))
% 
% 
% 	
% 
% History:
% Copyright Tommy Hosman, All Rights Reserved
%   2019   



% Force data to be a cell
if ~iscell(data)
    data = {data};
end


%% Handle mode and special cases
warnOnRemove = false; 
if ~exist('removeOutOfBounds','var')
    warnOnRemove = true; % Print if removeOutOfBounds was set implicitly
    removeOutOfBounds = 1;
end

if nargin < 4 || isempty(mode)
    mode = 'td';
end
if ~iscell(mode)
    mode = {mode};
end

if all(strcmpi(mode, 'dpca')) && length(data) == 2
    % Assume, data = {data, labels}
    % data mode output is 'etd', [nEvents x nTime x nDim]
    % label mode output is 'ed', [nEvents x nLblDim]
    mode = {'etd', 'ed'};
end



% Repeat mode for each input element if not specified
if isscalar(mode) && length(mode) < length(data)
    mode = repmat(mode,length(data),1);
end



%% Get epoch window
epoch = GetWindowAroundEvents(events, window);


if removeOutOfBounds
    epoch = HandleOutOfBoundEpoch(data, epoch, warnOnRemove);
end

%% Parse and reshape the inputs
for ii = 1:length(data)
    varargout{ii} = HandleDataElement(data{ii}, epoch, window, lower( mode{ii} ));
end
varargout{end+1} = epoch;


end



function data = HandleDataElement(data, epoch, window, mode)
    if isvector(data)
        data = data(:);
    end
    [nRows, nDim]    = size(data);
    [nTime, nEvents] = size(epoch);

    isDataEventBased   = nRows == nEvents;
    isOutputEventBased = strcmpi(mode,'ed');
    
    % Reshape to a default format of: [nTime x nEvents x nDim]
    if isDataEventBased
        if ~isOutputEventBased
            % If we are not outputting by event, then repeat the data
            % nEvents x 1 -> nTime x nEvents x nDim
            data = repmat( reshape(data, [1 size(data)]), nTime, 1, 1 );
        else
            % Data is nEvents x nDim and our output is nEvents x nDim
        end
    else
        
        
        
        try
            % Get epoch of data
            data = reshape(data(epoch(:),:), [size(epoch) nDim]);
            
            
        catch err
            % Catch the error to alert the user if due to epoch out of bounds
            me = max(epoch(:));
            mi = min(epoch(:));
            if me > nRows || mi<1
                fprintf(2, '\n\nERROR\nEvents + window is outside the bounds of data!\n');
            end
            rethrow(err)
        end
    end



    % Unless data is eventVector and mode is 'ed', data shape is expected to be
    % [nTime x nEvents x nDim]
    if ~(isDataEventBased && isOutputEventBased)
        expectedShape = [nTime, nEvents, nDim];
        dataShape = size(data,1:3);
        if ~isequal( dataShape, expectedShape)
            error('data in unexpected shape: [%s]. Expected: [%s]', num2str(dataShape), num2str(expectedShape));
        end
    end
    
    
    
    %% Reshape data to output mode
    
    switch mode
        case 'td' % time-event x dim
            data = reshape(data, nTime*nEvents, nDim);


        case {  'det' % defaultOrder(perms(1:3))
                'dte'
                'edt'
                'etd'
                'tde'
                'ted'}
            % Use ismember to map the default letters to mode
            defaultOrder = 'ted'; % default dim letter order
            [~,permInds] = ismember(mode,defaultOrder);
            data = permute(data, permInds);
            

        case 'ed' % Event x dim

            if isDataEventBased
                % Do nothing
            else
                % Choose which element to grab
                % Try to find the element at the onset of the event index 
                % (i.e., when window==0).
                % Otherwise take the first element.
                %
                % Possibly take the first element >= 0?
                eventInd = find(window==0, 1);
                if isempty(eventInd)
                    eventInd = 1;
                end
                data = squeeze(data(eventInd,:,:));
            end

        otherwise
            error('Unrecognized mode: %s', mode)
    end


end


function epoch = HandleOutOfBoundEpoch(data, epoch, warnOnRemove)
    %%
    [nTime, nEvents] = size(epoch);
    rmvEvents = false(1,nEvents);
    
    for ii = 1:length(data)
        d = data{ii};
        if isvector(d)
            d = d(:);
        end
        [nRows, nDim]    = size(d);
        isDataEventBased = nRows == nEvents;
        
        isOutOfBounds = epoch < 1;
        if ~isDataEventBased
            % Only check > nRows if data is step based
            isOutOfBounds = isOutOfBounds | epoch>nRows;
        end
        
        rmvEvents = rmvEvents | any(isOutOfBounds);
    end
    
    if any(rmvEvents) && warnOnRemove
        fprintf('Warning: Removing %d events + window that were of bounds.\n', sum(rmvEvents))
    end
    
    % Remove
    epoch(:, rmvEvents) = [];
end
