function h = InfoString( str, loc )
% Creates info annotation
% See SessionAnalysisInfo for standard analysis info string
%
%--------------------------------------------------------------------------
% History:
%   2021.06   Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------


%%
if nargin < 1
    str = '';
end
if nargin < 2
    loc = [0.005 0.7 0.15 0.3];
end
fontName = get(gcf,'defaultAxesFontName');

h = annotation('textbox', loc, 'string', str, 'edgecolor', 'none', 'FontName', fontName );

h.Interpreter = 'none'; % Default to none until we use tex

if ~nargout
    clear h
end
end