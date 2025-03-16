function h = YLabelString( str, loc, overwrite )
% h = YLabelString( str, loc, overwrite )
% 
% Creates a single ylabel on the left center of a figure.
%
% 
%--------------------------------------------------------------------------
% History:
%   2021.06   Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------


%%
tagStr = 'YLabelStr';
fh = gcf;
if nargin < 1
    str = '';
end
if nargin < 2
    x = GetXPos();
    loc = [x 0.5]; % ~ 0.05 X for tight_subplots
end
if nargin < 3
    overwrite = true;
end
if overwrite
    delete(findall(fh, 'Tag', tagStr));
end

fontSize = 20;
fontName = get(gcf,'defaultAxesFontName');
%%
% try; delete(ax); end
% try; delete(tbh); end

[X, Y] = GetXY(str, loc, fontSize);

% TextArrow is not always intuative...
% Can be helpful to add the headstyle and linestyle if debugging.
%
% Horizontal alignment to center gives different spacing in the X axis
% depending on the size of the text. HorizontalAlignment', 'left' gives
% consistent spacing, but the text is not centered per new line...
h = annotation(fh,'textarrow', X, Y, 'String',str, ...
    'FontSize',fontSize, 'FontWeight','bold', 'FontName', fontName,...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
    'HeadStyle', 'none', 'LineStyle', 'none', ...
    'TextRotation',90);

h.Tag = tagStr;
h.UserData.str = str;
h.UserData.loc = loc;
h.UserData.X = X;
h.UserData.Y = Y;
h.UserData.overwrite = overwrite;

%% On resize adjust the position again

fh.ResizeFcn = {@FigResizeCallback, h};
%%
if nargout == 0
    clear h
end
end


function x = GetXPos()
    %%
    ch = get(gcf, 'Children');
    minX = inf;
    for ii =1:length(ch)
        if strcmp(ch(ii).Type, 'axes')
            minX = min([ minX ch(ii).Position(1)]);
        end
    end
    x = minX - minX*0.3;
end


function [X, Y] = GetXY(str, loc, fs)
% Create axes and text object to get extent size (axes is required because
% text object in subplot will give incorrect values)
axPos = [0.1 0.1 .8 0.8];
ax = axes('Units', 'normalized', 'Position', axPos, 'OuterPosition', axPos,  'Visible', 'on', 'Interactions', []);
tbh = text(ax,0.5, 0.5, str, 'FontSize', fs, 'FontWeight', 'bold', 'Visible', 'on', 'Units', 'normalized');
hw = tbh.Extent(3:4);

try; delete(ax); end
try; delete(tbh); end
% try; delete(h); end


%% Create the y str label

X = [loc(1) 0];
Y = [loc(2) loc(2)];

% Adjust X, Y to have centered text based on the text extent
fudgeFactor = 0.01; % 0.01; % Found by calling with a single letter at 0.5 and shifting to be at 0.5;
Y = max(Y - hw(1)/2 + fudgeFactor, 0); % 
X(1) = max(X(1) - hw(2)/4, 0);
end


function FigResizeCallback(obj,event, h)
persistent inCallback
if ~isempty(inCallback)
    return
end
inCallback = true;

% Only adjust if defaults were used (position was not changed outside of
% this function)
if ~isvalid(h)
    obj.ResizeFcn = [];
    inCallback = [];
    return;
end

if isequal(h.X, h.UserData.X) && isequal(h.Y, h.UserData.Y)
    
%     YLabelString( h.UserData.str, h.UserData.loc, 0 )
    
    [X, Y] = GetXY(h.UserData.str, h.UserData.loc, h.FontSize);
    h.X = X;
    h.Y = Y;
    h.UserData.X = X;
    h.UserData.Y = Y;
    
    % Let everything else adjust before exiting
    pause(0.001); drawnow;
end


inCallback = [];


end