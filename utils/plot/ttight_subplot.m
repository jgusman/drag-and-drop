function [ha, pos, args] = ttight_subplot(nRow, nCol, varargin)
% ttight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = ttight_subplot(nRow, nCol, gap, marg_h, marg_w)
%
%   in:  nRows      number of axes in hight (vertical direction)
%        nCols      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 
%
%       mergeRC: {rowInds, colInds}
%           Cell array of row, col inds to merge.
%           returned handles still have the same nRows, nCols shape.
%           Can be a cell array of cell arrays, e.g.,
%           Ex1 tsp = ttight_subplot(4,4, 'mergeRC', {{1:2, 1}, {1:2, 3:4}}); % Merges left col and upper right;
%           Ex2 tsp = ttight_subplot(4,4, 'mergeRC', {1, 1:4}); % Merges top row. 
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% 
% TODO
%   Allow a set default argument that will set an appdata value with
%   updated defaults
% 
% Tommy Hosman 12.10.2021
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


%%
inputArgOrder = {'gap', 'marg_h', 'marg_w', 'proportional_h', 'proportional_w', 'mergeRC','xTickOffset','yTickOffset'};

dfArg.gap = 0.07; % 0.08; %0.05; 
dfArg.marg_h = [0.08 0.06]; %0.08; %0.06;
dfArg.marg_w = 0.06; %0.08; %0.06;
dfArg.proportional_h = []; % nRows that add up to 1
dfArg.proportional_w = []; % nCols that add up to 1
dfArg.mergeRC = {}; % Specify a subset of axes that are merged. ttight_subplot(2,2, 'mergeRC', {1, 1:2}). Merged axes all point to the same axes handle.
dfArg.xTickOffset = 1; % (matlab default 2) Tick label location: Negative is up, positive is down.
dfArg.yTickOffset = 0;  % (matlab default 2) Tick label location: Negative is right, positive is left.


%% Handle special case for arg hints (see functionSignatures.json)
if strcmp(nRow, 'keyvals')
    ha = inputArgOrder;
    return;
end


%% Normal case

args = HandleArgs(nRow, nCol, dfArg, inputArgOrder, varargin{:});

% Normalize
args.proportional_h = args.proportional_h / sum(args.proportional_h);
args.proportional_w = args.proportional_w / sum(args.proportional_w);


axh = (1-sum(args.marg_h)-(nRow-1)*args.gap(1)) .* args.proportional_h; 
axw = (1-sum(args.marg_w)-(nCol-1)*args.gap(2)) .* args.proportional_w;



% Plan position
% 
% Top to bottom
% Left to right
py = 1-args.marg_h(2); 
axPos = cell(nRow,nCol);
for ih = 1:nRow
    py = py-axh(ih);
    px = args.marg_w(1);
    for ix = 1:nCol
        axPos{ih, ix} = [px py axw(ix) axh(ih)];
        px = px+axw(ix)+args.gap(2);
    end
    py = py-args.gap(1);
end



%% Create axes
ha = gobjects(nRow,nCol);

% Merge any
if ~isempty(args.mergeRC)
    [mergeAx,mergeLoc] = MergeAxes(axPos, args);
    ha(mergeLoc) = mergeAx; % Copy ax handle to each merged axes
end


% Default
for ih = 1:nRow    
    for ix = 1:nCol
        if isa(ha(ih,ix), 'matlab.graphics.GraphicsPlaceholder')
            ha(ih,ix) = CreateAxesLocal(axPos{ih,ix},args);
        end
    end
end


if nargout > 1
    pos = get(ha,'Position');
end

% Set first axes to be focused
axes(ha(1));

% ha = ha(:);
% ha = reshape(ha,nCols,nRows)';


end

function ax = CreateAxesLocal(position, args)
    ax = axes('Units','normalized','Position',position);
    % Default gap offset value is 2
    ax.XRuler.TickLabelGapOffset = args.xTickOffset; %-4; % Tick label location: Negative is up, positive is down.
    ax.YRuler.TickLabelGapOffset = args.yTickOffset; %0; % Tick label location: Negative is right, positive is left.
end


function [mergeAx,mergeLoc] = MergeAxes(axPos, args)
    
    % Check if a cell of cells was passed
    if iscell(args.mergeRC{1}) && length(args.mergeRC{1}) == 2
        [mergeAx,mergeLoc] = deal([]);
        mergeRC = args.mergeRC;
        for ii = 1:length(mergeRC)
            args.mergeRC = mergeRC{ii};
            [tmpAx,tmpLoc] = MergeAxes(axPos, args);
            mergeAx  = cat(1, mergeAx, tmpAx(:));
            mergeLoc = cat(1, mergeLoc, tmpLoc(:));
        end
        return        
    end
    
    % Default case
    
    mR = args.mergeRC{1};
    mC = args.mergeRC{2};
    minmaxX = [inf -inf];
    minmaxY = [inf -inf];
    isMerge = false(size(axPos));
    for ri = 1:length(mR)
        for ci = 1:length(mC)
            isMerge(mR(ri),mC(ci)) = true;
            tmp = axPos{mR(ri),mC(ci)};
            tmpWH = tmp(1:2) + tmp(3:4);
            minmaxX(1) = min(minmaxX(1), tmp(1));
            minmaxX(2) = max(minmaxX(2), tmpWH(1));
            
            minmaxY(1) = min(minmaxY(1), tmp(2));
            minmaxY(2) = max(minmaxY(2), tmpWH(2));
        end
    end

mergePos = [minmaxX(1), minmaxY(1), diff(minmaxX), diff(minmaxY)];
ax = CreateAxesLocal(mergePos, args);
mergeAx  = repmat(ax,sum(isMerge(:)),1);
mergeLoc = find(isMerge);
end

function args = HandleArgs(nRows, nCols, dfArg, inputArgOrder, varargin)
    if ~isempty(varargin) && isstruct(varargin{1})
        opt.noMergeIfEmpty = true;
        args = MergeBintoA(dfArg, varargin{1}, opt);
    else
        args = dfArg;

        strVals = find(cellfun(@ischar, varargin));
        isAnyKeyVal = find(ismember(varargin(strVals), inputArgOrder));
        if ~isempty(isAnyKeyVal)
            for ii = 1:length(isAnyKeyVal)
                ind = strVals(isAnyKeyVal(ii));
                args.(varargin{ind}) = varargin{ind+1};
            end
        else        
            for ii = 1:length(inputArgOrder)
                if length(varargin) >= ii
                    % Input order
                    if ~isempty(varargin{ii})
                        args.(inputArgOrder{ii}) = varargin{ii};
                    end
                end
            end
        end
    end




if isscalar(args.gap)
    args.gap = [args.gap args.gap];
end
if isscalar(args.marg_w)
    args.marg_w = [args.marg_w args.marg_w];
end
if isscalar(args.marg_h)
    args.marg_h = [args.marg_h args.marg_h];
end
if isempty(args.proportional_h)
    args.proportional_h = ones(1,nRows)./nRows;
elseif length(args.proportional_h) == nRows-1
    args.proportional_h(end+1) = 1-sum(args.proportional_h);
end
if isempty(args.proportional_w)
    args.proportional_w = ones(1,nCols)./nCols;
elseif length(args.proportional_w) == nCols-1
    args.proportional_w(end+1) = 1-sum(args.proportional_w);
end

end