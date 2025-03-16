function SetTickLabels(varargin)
% SetTickLabels(xy,labels,tickVals,rot)
% SetTickLabels(ax,xy,labels,tickVals,rot)
% SetTickLabels(labels,tickVals,rot) % xy defaults to 'x'
% 
% Sets x or y tick labels
%
% 
% Inputs:
%   labels (required)
%       (cell str, categorical) The labels to use for the ticks.
%   ax (optional)
%       (axes handle) The axis to update.
%   xy (optional)
%       (char) 'x', 'y', 'xy', 'yx'. Specifies updating x, y or both ticks.
%   tickVals (optional)
%       Where to place the labels. Defaults to 1:length(labels).
%   rot (optional)
%       Rotation of tick labels in degrees.
%   
%--------------------------------------------------------------------------
% History:
%   2022.10      Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------


    [xy,labels,tickVals,ax, rot] = HandleInputArgs(varargin{:});
    
    for ii = 1:length(xy)
        switch lower(xy(ii))
            case 'x'

                ax.XTick = tickVals;
                ax.XTickLabel = labels;
                rot = ShouldRotateX(labels,rot);
                if ~isempty(rot)
                    ax.XTickLabelRotation = rot;
                end
            case 'y'
                ax.YTick = tickVals;
                ax.YTickLabel = labels;
                if ~isempty(rot)
                    ax.YTickLabelRotation = rot;
                end
        end
    end
    
end



function [xy,labels,tickVals,ax, rot] = HandleInputArgs(varargin)
    
    % ax
    isFirstAxes = ((isnumeric(varargin{1}) && isscalar(varargin{1}) && isgraphics(varargin{1}, 'axes')) ...
               || isa(varargin{1},'matlab.graphics.axis.AbstractAxes') || isa(varargin{1},'matlab.ui.control.UIAxes'));
    if isFirstAxes
        ax = varargin{1};
        varargin(1) = [];
    else
        ax = gca;
    end
    
    
    % xy
    if ischar(varargin{1})
        % First is xy
        xy = varargin{1};
        varargin(1) = [];
    else
        xy = 'x';
    end
    
    
    % labels
    labels = varargin{1};
    varargin(1) = [];
    
    
    % tickVals
    if isempty(varargin)
        tickVals = []; 
    else
        tickVals = varargin{1};
        varargin(1) = [];
    end
    if isempty(tickVals)
        tickVals = 1:length(labels);
    end
    
    % rot
    if isempty(varargin)
        rot = [];
    else
        rot = varargin{1};
        varargin(1) = [];
    end
end

function rot = ShouldRotateX(labels,rot)
    xRotateLabelThresh = 5;
    xRotateVal = 25;
    if iscategorical(labels)
        labels = cellstr(labels);
    end
    maxLen = max(cellfun(@length,labels));
    if length(labels) > xRotateLabelThresh && maxLen > 10 && isempty(rot) 
        rot = xRotateVal;
    end
end