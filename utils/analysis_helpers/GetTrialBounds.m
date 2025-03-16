function trialBounds = GetTrialBounds(in, type, relative)
% If relative, the first value is 1

[nR,nC] = size(in);
if nargin < 2
    if nC == 2
        type = 'startStops';
    else
        type = 'epoch';
    end
end

if nargin < 3
    relative = 1;
end


boundFun = @(x) cumsum( [1; diff( x,  [], 2 )+1 ] );
switch type
    case 'startStops'
        trialBounds = boundFun(in);
    case 'epoch'
        in = in([1 end],:)';
        trialBounds = boundFun(in);
    otherwise
        error('Unsupported %s', type)
end

if ~relative
    trialBounds = trialBounds + in(1) - 1;
end

end