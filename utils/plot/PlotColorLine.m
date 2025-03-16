function h = PlotColorLine(varargin)
% Plots a 2D or 3D line with interpolated color using a patch object.
% 
% h = PlotColorLine(x,y,c,varargin)
% h = PlotColorLine(x,y,z,c,varargin)
% 
% h = PlotColorLine(xy,c,varargin)
% h = PlotColorLine(xyz,c,varargin)
% 
% h = PlotColorLine(xy,varargin)     % Generates c as linear vect across line 
% h = PlotColorLine(xyz,varargin)    % Generates c as linear vect across line 
% 
% h = PlotColorLine(ax,...,varargin) % First arg can always be an axes handle
% h = PlotColorLine(); % Plots demo - see examples
% 
% 
% https://www.mathworks.com/help/matlab/ref/patch.html
% https://www.mathworks.com/help/matlab/creating_plots/add-transparency-to-graphics-objects.html
% 
% 
% Example1
% x = linspace(-pi,pi,100);
% randXY = cat(1,x,cos(x))';
% PlotColorLine(randXY)
% 
% Example2
% c = cos(x);
% PlotColorLine(randXY, c, 'linewidth', 3)
% 
% Example3
% x = linspace(-pi,pi,100);
% randXYZ = cat(1,x,cos(x),sin(x))';
% PlotColorLine(randXYZ)
% view(3)
% 
% See also patch
% -------------------------------------------------------------------------

% No args, demo (examples)
if isempty(varargin)
    tmpH=RunDemo(); 
    if nargout > 0, h=tmpH; end
    return;
end
    

    [ax, xyz, c, patchArgs] = GetInputs(varargin{:});
    
    tmpH = patch(ax, xyz{:},  c, ...
                 'EdgeColor','interp',...
                 'MarkerFaceColor','flat', ...
                 patchArgs{:});
    %%
    if nargout > 0
        h = tmpH;
    end
end




function [ax,xyz,c,vOut] = GetInputs(varargin)
% h = PlotColorLine(x,y,c,varargin)
% h = PlotColorLine(xy,c,varargin)
% h = PlotColorLine(xy,varargin)
% h = PlotColorLine(xyz,...,varargin)
% h = PlotColorLine(ax,...,varargin)


% An axes
if ~isempty(varargin) && isscalar(varargin{1}) && ishghandle(varargin{1}, 'axes')
    ax = varargin{1};
    varargin(1) = [];
else
    ax = gca;
end


% Find how (x), y (z) were passed
if ~isempty(varargin)
    
    nArgs = length(varargin);
    nVectInputs = 0;
    
    if nArgs > 2
        % Count how many numeric vectors there are.
        % If 2, assume x,y
        % If 3, assume x,y,c
        % If 4, assume x,y,z,c
        
        for ii = 1:nArgs
            if isnumeric(varargin{ii}) && isvector(varargin{ii})
                % Check if third argument are all integers
                % Assume if a vector of integers, that this is a color
                if ii==3
                    isAllIntegers = isequaln(varargin{ii},round( varargin{ii} ));
                    if isAllIntegers
                        % Color!
                        break;
                    end
                end
                
                % This is an input to be plotted
                nVectInputs = nVectInputs + 1;
            else
                % Exit on first non numeric vector
                break
            end
        end
    end
    
    if ismember(nVectInputs, 1:3)
        % x, y, (z) are input 1, 2, (3)
        xyz = varargin(1:nVectInputs);
        varargin(1:nVectInputs) = [];
        
    elseif ismember(nVectInputs, 4)
        % x,y,z are inputs 1-3, color is a vector
        xyz = varargin(1:3);
        varargin(1:3) = [];

    elseif ~isvector(varargin{1})
        % xy (xyz) is a 2 (3) dim matrix
        
        nDim = size(varargin{1},2);
        if ismember(nDim,[2,3])
            for ii = 1:nDim
                xyz{ii} = varargin{1}(:,ii);
            end
            varargin(1) = [];
        else
            error('Unknown condition 1');
        end    
    else
        error('Unknown condition 2');
    end
else
    error('Not enough inputs. See the help for input methods')
end

xyzLen = cellfun(@length, xyz);
for ii = 1:length(xyz)
    if xyzLen(ii) == 1
        xyz{ii} = repmat( xyz{ii}, max(xyzLen), 1);
    end
end

% Color
c = 1:numel(xyz{1}); % Default

if ~isempty(varargin)
    if  (isnumeric(varargin{1}) && isvector(varargin{1})) || ... % [1 2 3 4]
        (isnumeric(varargin{1}) && size(varargin{1},2) == 3) % [r,g,b]
        
        % color as a char is not valid 

        c = varargin{1};
        varargin(1) = [];
    end
end


% Misc clean up
if isvector(c)
    c = c(:);
end

% Prevent the closing of the patch object on itself
for ii = 1:length(xyz)
    xyz{ii}(end+1) = nan;
end

c(end+1,:) = c(end,:);

if ~isvector(c) && isnumeric(c) && size(c,ndims(c)) == 3
    % a color matrix of shape [N x 3] needs to be a [1 x N x 3]
    c = reshape(c, [1 size(c)]);
end

% Additional line arguments
vOut = varargin;
end

function h = RunDemo()

figure
subplot(311)
% Example1
x = linspace(-pi,pi,100);
xy = cat(1,x,cos(x))';
h(1) = PlotColorLine(xy);
title(sprintf('{\\fontsize{14}PlotColorLine Demo}\nExample 1'))




subplot(312)
% Example2
c = cos(x);
h(2) = PlotColorLine(xy, c, 'linewidth', 3);
title('Example 2')


subplot(313)
% Example3
x = linspace(-pi,pi,100);
xyz = cat(1,x,cos(x),sin(x))';
h(3) = PlotColorLine(xyz);
view(3); % 3D
title('Example 3')
end

