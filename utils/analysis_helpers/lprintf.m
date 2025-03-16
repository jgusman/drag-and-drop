function nOut = lprintf(varargin)
% nOut = lprintf(varargin)
% Used for printing output that overwrites the previous output, e.g.
% printing percent complete in a loop.
% 
% Use like fprintf, but call without any inputs before the for loop to
% initialize (set 'prevPrintLen' to 0).
% 
% This is used because carriage return creates a new line in matlab.
% 
% Idea to concat chars with 1 print came from
% https://undocumentedmatlab.com/articles/command-window-text-manipulation
% 
% Example:
%     N = 10;
%     lprintf(); % Init
%     for ii = 1:N
%         lprintf('Loop %d out of %d\n', ii, N)
%         pause(0.1)
%     end
%--------------------------------------------------------------------------
% History:
%   2021.03      Copyright Tommy Hosman, All Rights Reserved
%--------------------------------------------------------------------------

    persistent prevPrintLen
    
    
    % Clear prevPrint
    if isempty(varargin)
        prevPrintLen = 0;
        return;
    end

    % If FID is passed, grab that
    if isnumeric(varargin{1})
        fid = varargin{1};
        varargin(1) = [];
    else
        fid = 1;
    end
    
    % Now create the backpace str and print str
    removePrev   = repmat(sprintf('\b'), 1, prevPrintLen);
    printStr     = sprintf(varargin{:});
    prevPrintLen = length(printStr);
    
    % Print
    fprintf(fid,[removePrev printStr]);
    
    % Output
    if nargout == 1
        nOut = prevPrintLen;
    end
end