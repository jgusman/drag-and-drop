function lr = LeftRightPlot(fh)
%     LeftRightPlot(fh)
% 
% plotFun expects the figure handles as the input
% plotFun(fh)
% 
% fh = figure;
% lr.ind = 1;
% lr.minMax = [1 length(outSet)];
% lr.plotFun = @PlotAIFromFigUserData;
% lr.wrapBounds = false;
%
% fh.UserData.lr = lr;
% fh.UserData.data = data;
% LeftRightPlot(fh)


% Default

lr = GetDefaultLR();

if nargin == 1
    % Make sure defaults are in lr
    if isfield(fh.UserData, 'lr')
        lr = MergeBintoA(lr,fh.UserData.lr);
    end
    fh.UserData.lr = lr;
    
    if isfield(fh.UserData, 'ud')
        fh.UserData.ud = MergeBintoA(GetDefaultLR(),fh.UserData.ud);
    end
    
    % Set up key press
    UseMyKeypress(fh, @UpdateLR_Callback)
end

%% Return default lr struct
if ~nargout
    clear lr
end
    
end

function lr = GetDefaultLR()
    lr.ind = 1;
    lr.minMax = [1 1];
    lr.plotFun = [];
    lr.wrapBounds = true;
    lr.inc = 1; % Increment amount
    lr.dec = 1; % Decrement amount
end

function UpdateLR_Callback(obj,evnt)
% Expects lr struct in userdata
% lr.ind
% lr.minMax
% lr.plotFun
    if isempty(obj.UserData)
        return
    end
    if ~isfield(obj.UserData, 'lr')
        error('This callback requires ''lr'' struct to be in UserData')
    end

    plotFun = obj.UserData.lr.plotFun;
    if ischar(plotFun)
        plotFun = str2func(plotFun);
    end

    switch evnt.Key
        case 'rightarrow'
            IncObj(obj,'lr',plotFun)
            
        case 'leftarrow'
            DecObj(obj,'lr',plotFun)
            
        case 'uparrow'
            IncObj(obj,'ud',plotFun)
        case 'downarrow'
            DecObj(obj,'ud',plotFun)
    end
    
end

function IncObj(obj,fn,plotFun)

    if isfield(obj.UserData,fn)
        newInd = obj.UserData.(fn).ind+obj.UserData.(fn).inc;
            
        if newInd > obj.UserData.(fn).minMax(2)
            if obj.UserData.(fn).wrapBounds
                % Wrap to start
                newInd = obj.UserData.(fn).minMax(1);
            else
                % Stay to end
                newInd = obj.UserData.(fn).minMax(2);
            end
        end
        obj.UserData.(fn).ind = newInd;
        plotFun(obj)
    end
end

function DecObj(obj,fn,plotFun)
    if isfield(obj.UserData,fn)
        newInd = obj.UserData.(fn).ind-obj.UserData.(fn).dec;
            
        if newInd < obj.UserData.(fn).minMax(1)
            if obj.UserData.(fn).wrapBounds
                % Wrap to end
                newInd = obj.UserData.(fn).minMax(2);
            else
                % Stay at start
                newInd = obj.UserData.(fn).minMax(1);
            end
        end
        obj.UserData.(fn).ind = newInd;
        plotFun(obj)
    end
end
