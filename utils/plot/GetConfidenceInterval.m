function [u,ci] = GetConfidenceInterval(data,method)
% [u,ci] = GetConfidenceInterval(data,method)
    if nargin < 2
        method = 'normfit';
    end
    
    switch method
        case 'ci'
            nBoot = 1000;
            if isvector(data)
                data = data(:);
            end
            % Follow normfit, iterate of dim 2
            for ii = 1:size(data,2)
                ci(:,ii) = bootci(nBoot, {@nanmean, data(:,ii)});
            end
            u  = nanmean(data,1);
            ci = [u       - ci(1,:); ...
                  ci(2,:) - u];
             %%
            
        case 'normfit'
            if any(isnan(data(:)))
                d = data(~isnan(data));
            else
                d = data;
            end
            
            [u,~,ci] = normfit(d);
            if isscalar(d)
                ci = nan;
            else
                ci = u - ci(1,:);
            end
        otherwise
            error('Unsupported method: %s', method);
    end
    
    
    
end