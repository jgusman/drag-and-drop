function blkStr = BlocksToString(blocks,separator)
if nargin < 2
    separator = '-';
end

    blocks(isnan(blocks)) = [];
    if isempty(blocks)
        blkStr = '';
    elseif isscalar(blocks)
        blkStr = sprintf('%d ', blocks);
    elseif all( diff(blocks) == 1 )
        blkStr = sprintf('%d%s%d', min(blocks),separator, max(blocks));
    else
        % Handle single and consecutive blocks
        blocks = sort(blocks(:)'); % Make it a row
        
        blkDiff = diff(blocks);
        
        [~, consecutiveBlks] = FindConsecutiveOnes( blkDiff == 1,1);
        if ~isempty(consecutiveBlks)
            consecutiveBlks(:,2) = consecutiveBlks(:,2) + 1;
        end
        singleBlks = setdiff( blocks,  blocks( RowColon(consecutiveBlks)  ));
        singleBlks = repmat( singleBlks(:), 1, 2 );
        
        allBlocks = cat(1, singleBlks, blocks(consecutiveBlks));

        
        [~, si] = sort(allBlocks(:,1));
        allBlocks = allBlocks(si,:);
        
        blkStr = '';
        for ii = 1:size(allBlocks,1)
            if diff(allBlocks(ii,:)) == 0
                % Single block
                blkStr = sprintf('%s %d, ', blkStr, allBlocks(ii,1) );
            else
                blkStr = sprintf('%s %d%s%d, ', blkStr, allBlocks(ii,1),separator,allBlocks(ii,2) );
            end
        end
        
        blkStr(end-1:end) = [];
        blkStr = strtrim(blkStr);
    end
end