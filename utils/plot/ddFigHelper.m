classdef ddFigHelper
% ddFigHelper
% 
% Static class to help export paper quality figures.
% 
% There seems to be no perfect method to export high quality figures.
% Note that I've had the best luck with .svg output for vector files
% (however this still requires re-scaling with post-process program, e.g.
% illustrator).
% 
% Confirm your figure size and font size with external program (see demo)! 
% 
% Example
%   ddFigHelper.Demo()
% 
% Example2
%     ddFigHelper.SetSaveDir(fullfile(pwd,'FigEx2'))
%     figWidthHeight = [3 3]; % inches
%     ddFigHelper.CreateFigure(100,figWidthHeight); % Or use ddFigHelper.UpdateFigDims to set figure width/height
%     clf
%     plot(1:10)
%     text(0,0,'Text font should match defaults')
%     ddFigHelper.AddLetter('A')
%     ddFigHelper.SaveFigure('example2')
% 
% Methods:
% CreateFigure, SaveFigure, AddLetter
% 
% 
% Methods to update fig helper properties
%   figParams = ddFigHelper.GetAppData()
%   ddFigHelper.SetSaveDir( saveDir )
%   ddFigHelper.SetParams( figParams_new )
% 

    properties (Constant)
        % Can modify params here or in app data (see ddFigHelper.Demo())

        fig = struct('color',   [1 1 1], ...
                     'fontSize', 9,...
                     'fontName', 'Arial', ...
                     'units', 'inches');
        

        % Subpanel letter attributes         
        letters = struct('FontSize', 15, ...
                         'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'left', ...
                         'VerticalAlignment', 'bottom');
           
        
        % Must be specified to save
        % Can be specified with a hardcode here
        % or, by calling ddFigHelper.SetSaveDir('save dir path'), 
        % up calling ddFigHelper.SetParams( )
        saveDir = '';
        
        % See SaveFigure() for supported formats
        saveFormats = {'png','svg'};
        
        
        % Specific formats can be placed inside subdirectories
        saveFormatSubdirs = struct('svg','vectors');
        
               
        
        % Due to a matlab bug with svg output. 
        % Show a text message on the image that states the image should be
        % scaled to get exact font/figure specs.
        %
        % Message states expected figure size and reference font
        % see ddFigHelper.CreateScaleMessage()
        % List save formats to display this message.
        
        showScaleMessageOnSave = {'svg'}; 
        
        
        logSubdir = 'statistics'; % subdir of saveDir for Log
    end
    
    
    % Private app data name
    properties (Constant, Access=private)
        appdataName = 'ddFigHelper';
    end
    
    
    
    
    %% Methods, Create, Save, Add letter
    
    methods (Static)
        function fh = CreateFigure(figNum, figWH)
            % fh = CreateFigure(figNum, figWH)
            % Creates figure with default fonts, axes
            if nargin < 1 || isempty(figNum)
                fh = figure();
            else
                fh = figure(figNum);
            end
            if nargin < 2
                figWH = []; % figure width, height In units as specified in fig.units
            end
            
            figParams = ddFigHelper.GetAppData();
            sFig = figParams.fig;
            
            % If we use set str names in ddFigHelper.fig, it can be
            % arbitrarily set without updating below (like letter).
            %
            % But, it is not fun to reference DefaultTextFontName elsewhere
            set(fh, ...
                'Color', sFig.color, ...
                'DefaultAxesFontName', sFig.fontName, ...
                'DefaultAxesFontSize', sFig.fontSize, ...
                'DefaultTextFontName', sFig.fontName, ...
                'DefaultTextFontSize', sFig.fontSize, ...
                'DefaultAxesFontSizeMode', 'manual', ...
                'Units', sFig.units, ...
                'DefaultAxesNextPlot','add'...  % essentially, permenant "hold on"
                );

            if ~isempty(figWH)
                ddFigHelper.UpdateFigDims(figWH,fh);
            end
            
            %% output
            if ~nargout
                clear fh
            end
        end
        
        function SaveFigure( fh, savename )
            % SaveFigure( fh, savename )
            % Saves figure with formats from ddFigHelper.saveFormats
            %
            % Note ddFigHelper.saveDir must be set (not empty) or savename
            % must be a full path to save.
            %
            % Do we want to assume that a savename with folders should
            % ignore saveDir? 
            % This will error if both savename and saveDir are filled out.
            
            if nargin < 2 && ~isa(fh,'matlab.ui.Figure')
                savename = fh;
                fh = gcf;
            end
            
            if strcmp(savename,'clip')
                saveFileType = '-dbitmap';
                orgRender = fh.Renderer;
                fh.Renderer = 'Painters';
                print( fh, '-clipboard', saveFileType )
                fh.Renderer = orgRender;
                return
            end
                
                
            % Get params
            figParams = ddFigHelper.GetAppData();
            saveDir = figParams.saveDir;
            
            
            % If savename is a full path, default to using that?
            ind = find(ismember(savename,filesep),1,'last'); 
            fileDir = savename(1:ind); % Can handle some paths better than filesep
            if ~isempty(fileDir)
                % And require that the current saveDir is empty?
                if isempty( figParams.saveDir )
                    ddFigHelper.SetSaveDir( fileDir );
                    saveDir = fileDir;
                end
            end
            
            
            % Check if saveDir is not empty
            if ~isempty(saveDir)

                
                %%
                
                fh.Name = sprintf('%s', savename);
                saveFormats = figParams.saveFormats;
                if ~iscell(saveFormats)
                    saveFormats = {saveFormats};
                end
                
                
                
            
                % Loop through to save each requested format
                for ii = 1:length(saveFormats)
                    % Creates a message to warn user about necessary scaling
                    % Will only create/display message if saving format
                    % is in ddFigHelper.showScaleMessageOnSave
                    ddFigHelper.CreateScaleMessage(fh,saveFormats{ii});
                    
                    % Create save path
                    savePath = ddFigHelper.CreateSaveSubdir(saveFormats{ii}, savename);
                    
                    switch saveFormats{ii}
                        case 'svg'
                            % If print does not work for some reason, try
                            % https://github.com/jschwizer99/plot2svg

                            savePath = [savePath '.svg'];
                            
                            saveFileType = '-dsvg';
                            orgRender = fh.Renderer;
                            fh.Renderer = 'Painters';
                            print( fh, savePath, saveFileType )
                            fh.Renderer = orgRender;
                            

                        case 'png'
                            dotsPerInch = '-r300';
                            savePath = [savePath '.png'];
                            saveFileType = '-dpng';
                            orgRender = fh.Renderer;
                            fh.Renderer = 'Painters';
                            print( fh, savePath, saveFileType, dotsPerInch )
                            fh.Renderer = orgRender;
                            
                        case 'pdf'
                            % NOTE font name may not perfectly be mapped to
                            % all elements for pdf
                            % see export_fig for non transparent pdf help
                            % use svg if you need transparency
                            
                            savePath = [savePath '.pdf'];
                            saveFileType = '-dpdf';
                            orgRender = fh.Renderer;
                            fh.Renderer = 'Painters';
                            print( fh, savePath, saveFileType )
                            fh.Renderer = orgRender;

                            
                            
                        otherwise
                            error('Unknown save format %s. Please add to ddFigHelper.SaveFigure()', saveFormats{ii})
                    end
                    
                    
                end
                
                % Clear scale message if it was present
                ddFigHelper.CreateScaleMessage(fh,'delete-scale-message');
                ddFigHelper.SetSaveDir( figParams.saveDir );
                
                fprintf('Saved Fig: %s\n',savename)
                
            end
        end
        
        function UpdateFigDims( figWH, fh )
            % UpdateFigDims( figWH, fh )
            %
            % Updates figure to have width and height of figWH
            % figWH = [width height] in units ddFigHelper.fig.units
            if nargin < 2
                fh = gcf;
            end


            %%
            figParams = ddFigHelper.GetAppData();
            figUnitsOrg = fh.Units;
            
            
            % I think this works best to not scale
            % and scale after saving the fig. The reason for this is that
            % font size uses pt. So, if we rescale the figure, font sizes
            % will not be appropriate (possibly due to a matlab svg bug).
            
            % scalePtPx = ddFigHelper.GetFigPixlePointRatio(fh);

            fh.Units = figParams.fig.units;
            fh.PaperUnits = figParams.fig.units;
            fh.Position(3:4) = figWH; %.*scalePtPx;
            fh.PaperPosition(3:4) = figWH; %.*scalePtPx;

            %%
            
            
            % Recenter figure if off the desktop
            % Otherwise, top right of figure (with close buttons) may be
            % off screen.
            ddFigHelper.AdjustFig(fh);
            
            fh.Units = figUnitsOrg;
            
        end
        
        function axArgs = GetDefaultAxArgs()
            % see ttight_subplot.m for how these are used
            axArgs.gap = 0.1; % 0.08; %0.05; 
            axArgs.marg_h = [0.1 0.08]; %0.08; %0.06;
            axArgs.marg_w = 0.06; %0.08; %0.06;
            axArgs.xTickOffset = 0; % (matlab default 2) Tick label location: Negative is up, positive is down.
            axArgs.yTickOffset = 1; % (matlab default 2) Tick label location: Negative is right, positive is left.
            axArgs.mergeRC = {}; 
            axArgs.proportional_h = [];
            axArgs.proportional_w = []; %[0.1 0.9]; % nRows that add up to 1
        end
        
        
        
        
        
        function letter_h = AddLetter(ax,letter,xOffset,yOffset)
            % an = AddLetter(ax,letter,xOffset,yOffset)
            %
            % Creates letter annotation. Places at top left of the passed
            % axes (or current axes if not passed).
            %
            % letter may be the character, e.g., 'a' or 'A'. 
            % if letter is a number, it is interpreted as 1 indexed into
            % the alphabet. For example, AddLetter(ax,1) places an 'A'
            %
            % x,y offsets are in normalized figure units
            % Maybe this should be normalized axes units?
            % Note that yOffset is from the top
            if ~exist('letter', 'var')
                letter = [];
            end
            if ~exist('xOffset','var') || isempty(xOffset)
                xOffset = -0.06; % Normalized from start of ax
            end
            if ~exist('yOffset','var') || isempty(yOffset)
                yOffset = -0.01; % Normalized from top of ax
            end
            if isempty(ax)
                ax = gca;
            end
            if ~isa(ax,'matlab.graphics.axis.Axes')
                if ~isempty(letter)
                    xOffset = letter;
                end
                letter = ax;
                ax = gca;
            end
            
            
            if isnumeric(letter) && isscalar(letter) && letter < 65
                numOffset = double('A'); % Use for capital letter index
                % numOffset = double('a'); % Use for lower case letter index
                letter = char( (numOffset-1) + letter ); % 64 is 65-1, 65 is 'A'
            end
            
            
            axU = ax.Units;
            ax.Units = 'normalized';
            pos = ax.Position;
            ax.Units = axU;
            
            figParams = ddFigHelper.GetAppData();
            
            
            letWH = [0.05 0.05];
            letXY = [pos(1)+xOffset  pos(2)+pos(4)+yOffset];
                        
            letter_h = annotation('textbox',[letXY letWH],'String',letter,'LineStyle','none', 'FontName', figParams.fig.fontName);
            
            % Set arguments
            % Equivalent to:
            % letter_h.FontSize = 12;
            letterArgs = figParams.letters;
            flatArgs = cat(2,fieldnames(letterArgs), struct2cell(letterArgs))';
            flatArgs = flatArgs(:);
            set(letter_h,flatArgs{:});
            
            
            % Delete old letter if it existed
            fh = ax.Parent;
            if isfield(fh.UserData, (letter))
                delete( fh.UserData.(letter) )
            end
            fh.UserData.(letter) = letter_h;
            
            %%
            if ~nargout
                clear letter_h
            end
            
        end
        
        function LogPrints(logfile)
            % LogPrints(logfile)
            %
            % Usage
            % logfile = 'log_figure3.txt';
            % ddFigHelper.LogPrints(logfile)
            % some calcs and prints...
            % 
            % % Three ways to turn off log / diary
            % diary off; % matlab command
            % ddFigHelper.LogPrints(); % no input defaults to diary off
            % ddFigHelper.LogPrints('off'); % 'off' keyword turns diary off
            
            %% Finished logging / Turn off diary
            if nargin < 1
                diary off
                return
            end
            
            if strcmpi(logfile,'off')
                diary off
                return
            end
            
            
            %% Path Setup
            figParams = ddFigHelper.GetAppData();
            saveDir = figParams.saveDir;
            
            % Add .txt if no extension
            % Find extension. Better than fileparts
            ind = find(ismember(logfile,'.'),1,'last');
            extStr = logfile(ind:end);
            if isempty(extStr)
                logfile = [logfile '.txt'];
            end
            
            
            logPath = fullfile(saveDir,figParams.logSubdir,logfile);
            
            
            % If savename is a full path, default to using that?
            ind = find(ismember(logPath,filesep),1,'last'); 
            logDir = logPath(1:ind); % Can handle some paths better than filesep
            if ~isempty(logDir)
                if ~isfolder(logDir)
                    mkdir(logDir)
                end
            end
            
            
            %% Create log
            if ~isempty(logDir)
                
                % Delete existing log files (do not append)
                % But maybe we want to append?
                if isfile(logPath)
                    delete(logPath); 
                end
                
                % Close any open diary file
                diary off
                % Start log
                diary(logPath)
            end
            
        end
        
        
        
    end
    
    
    %% Method Demo
    methods (Static)
        function Demo()
            % Create a figure and plot difficult to export data
            % Make it 6" x 3"
            % Make default fig font 8
            % 
            
            
            %% Setup defaults (can also change const props)
            
            saveDir = fullfile(pwd,'DemoFigs');
            fprintf('Saving demo files in current matlab path at\n%s\n', saveDir)
            
            
            originalParams = ddFigHelper.GetAppData();
            
            ddFigHelper.ResetParams();
            figParams = ddFigHelper.GetAppData();
            
            
            figParams.saveDir = saveDir;
            figParams.fig.units = 'inches';
            figParams.fig.fontSize = 8;
            figParams.fig.color = [1 1 1];
            figParams.fig.fontName = 'Arial';
            % What formats to save
            figParams.saveFormats = {'png','svg','pdf'};
            
            % Create subdirs for specific output file type (.svg maps to
            % save format svg above)
            figParams.saveFormatSubdirs.svg = 'svg-files';
            
            % Update params/settings
            ddFigHelper.SetParams( figParams )
            
            
            %% Create figure
            figWH = [6 3];
            figNum = 192927;
            

            fh = ddFigHelper.CreateFigure(figNum, figWH);
            clf(fh);
            
            % Create 3 axes. First axes twice the width as the others
            axs(1) = subplot(1,4,1:2);
            axs(2) = subplot(1,4,3);
            axs(3) = subplot(1,4,4);

            % Shift first, last axes a bit
            axs(1).Position(1) = 0.07;
            axs(3).Position(1) = axs(3).Position(1) + 0.05;
            
            %% Axes 1
            % Transparent 3D surface
            % With text, and line plot
            
            ax = axs(1);
            axes(ax)
            cla(ax)
            hold on
            
            [X,Y,Z] = peaks(25);
            h = surf(X,Y,Z);
            h.FaceAlpha = 0.5;
            
            x = linspace(-5,5,40);
            y = sin(x);
            plot3(x,y,x,'r','LineWidth', 3)
            
            demoStr = sprintf('Demo \nfigure size: [%.1f x %.1f] %s', ...
                figWH(1),figWH(2), figParams.fig.units);
            text(0,6,6,demoStr, 'VerticalAlignment','bottom')
            
            fontStr = sprintf('font name: %s\nfont size: %.1f',  ...
                figParams.fig.fontName, figParams.fig.fontSize);
            th = text(0,6,6, fontStr, 'VerticalAlignment','top',  ...
                                      'Color', [0.9624, 0.4456, 0.2654], 'FontWeight', 'bold');
            
            
            xlabel('x label')
            ylabel('y label')
            zlabel('z label')
            title('Axes 1')
            view(3)
            
            % Plot letter with bigger x offset
            letterOffsetX = -0.07;
            ddFigHelper.AddLetter('A',letterOffsetX)
            
            
            %% Axes  2
            % Bar plot with alpha
            
            ax = axs(2);
            axes(ax)
            cla(ax)
            hold on
            bar(1:4,'FaceAlpha',0.5)
            xlabel('Bar label')
            ylabel('Bar value (units)')
            title('Axes 2')
            ddFigHelper.AddLetter('B')
            
            %% Axes  3
            % Line plot with markers and font info text
            
            axisIndex = 3;
            ax = axs(axisIndex);
            axes(ax)
            cla(ax)
            hold on
            plot(x,y,'-p')
            xlabel('x label')
            ylabel('y label')
            title('Axes 3')
            
            
            
            % Can add letter with index (1='A', or 3='C')
            ddFigHelper.AddLetter(axisIndex)
            
            
            %% Save
            ddFigHelper.SaveFigure('Demo fig')
            
            % Note, if your pixel to point ratio is not 1, you will need to
            % scale the final output to match the desired fig/font sizes.
            %
            % This may just be a matlab svg bug.
            %
            % See:
            % ddFigHelper.GetFigPixlePointRatio()
            %   Shows the ratio
            % ddFigHelper.showScaleMessageOnSave
            %   Plots warning text on output image that states how much to
            %   scale image.
            
            %% Revert to original fig settings
            ddFigHelper.SetParams( originalParams );
        end
    end
    
    
    %% Get Set methods
    methods (Static)
        function figParams = GetAppData()
            figParams = getappdata(0,'ddFigHelper');
            if isempty(figParams)
                ddFigHelper.ResetParams();
                figParams = getappdata(0,'ddFigHelper');
            end
        end
        function ResetParams()
            warning off;
            figParams = struct(ddFigHelper);
            warning on;
            setappdata(0,ddFigHelper.appdataName,figParams);
        end
        
        function SetSaveDir(saveDir)
            figParams = ddFigHelper.GetAppData();
            figParams.saveDir = saveDir;
            setappdata(0,ddFigHelper.appdataName,figParams);
        end
        
        function SetParams( figParams_new )
            % SetParams( figParams_new )
            %
            % sets fields that are in figParams_new
            % Does a simple, "flat" set of figParams
            %
            % Example
            %  figParams = ddFigHelper.GetAppData();
            %  figParams.fig.fontName = 'Times New Roman';
            %  ddFigHelper.SetParams( figParams )
            
            
            figParams = ddFigHelper.GetAppData();
            
            fn = fieldnames(figParams_new);
            for ii = 1:length(fn)
                figParams.(fn{ii}) = figParams_new.(fn{ii});
            end
            
            setappdata(0,ddFigHelper.appdataName,figParams);
        end
        
        
        function pxPtRatio = GetFigPixlePointRatio(fh)
            if nargin < 1
                fh = gcf;
            end
            
            figUnitsOrg = fh.Units;
            
            fh.Units = 'points';
            posPt = fh.Position(3:4);

            fh.Units = 'pixels';
            posPx = fh.Position(3:4);
            pxPtRatio = posPt./posPx;
            
            fh.Units = figUnitsOrg;
        end
    end
    
    
    
    
    
    %% Private helper methods
    
    methods (Static, Access = private )
        function savePath = CreateSaveSubdir(saveFormat, savename)
            figParams = ddFigHelper.GetAppData();
            savePath = figParams.saveDir;
            
            % format specific subfolder?
            if isfield(figParams.saveFormatSubdirs, saveFormat)
                savePath = fullfile(savePath, figParams.saveFormatSubdirs.(saveFormat));
            end
            
            % Create dir if does not exist
            if ~isfolder(savePath)
                mkdir(savePath)
            end
            
            % append file name to save path
            savePath = fullfile(savePath, savename);
        end
        
        
        
        
        function AdjustFig(fh)
            fhU = fh.Units;
            fh.Units = 'pixels';
            
            figEdge = fh.Position(1:2) + fh.Position(3:4);
            monitorEdge = ddFigHelper.GetMonitorEdge(fh);            
            
            
            offScreenPix = monitorEdge - figEdge;
%             adjustXY = offScreenPix .* (offScreenPix < 0); 
            adjustXY = [0 0];
            for ii = 1:length(offScreenPix)
                if offScreenPix(ii) < 0
                    adjustXY(ii) = offScreenPix(ii);
                end
            end
            fh.Position(1:2) = fh.Position(1:2) + adjustXY;
            fh.Units = fhU;
        end
        
        function monitorEdge = GetMonitorEdge(fh)
            gr = groot;
            monitorStarts = gr.MonitorPositions(:,1:2); 
            monitorEdges = gr.MonitorPositions(:,1:2) + (gr.MonitorPositions(:,3:4)-1);
            monitorCenters = monitorStarts + (monitorEdges+1)./2;
            
            nMonitors = size(gr.MonitorPositions,1);
            
            figCenter = fh.Position(1:2) + fh.Position(3:4)./2;
            figCenters = repmat(figCenter,nMonitors, 1);
            
            isOverlap = [figCenters >= monitorStarts  figCenters <= monitorEdges];
            
            % Which monitor is the figure center on?
            monitorIndex = find(all(isOverlap,2));
            
            if isempty(monitorIndex)
                % Find closest to monitor center
                absDistFromFigCenter = abs(sum(monitorCenters - figCenter,2));
                [~,monitorIndex] = min(absDistFromFigCenter);
            end
            monitorEdge = monitorEdges(monitorIndex,:);
        end
        
        
        function CreateScaleMessage(fh, saveFormat)
            
                scaleMsgTag = 'PixPtWarning';
                delete(findall(fh,'Tag',scaleMsgTag))
                
                figParams = ddFigHelper.GetAppData();
                
                shouldShowRatio = any(strcmpi(figParams.showScaleMessageOnSave, saveFormat));
                if ~shouldShowRatio
                    return
                end
                        
                
            
                pxPtRatio = ddFigHelper.GetFigPixlePointRatio(fh);
                mismatchedRatio = any(pxPtRatio ~= 1);
                    
                if mismatchedRatio
                    mismatchedFontSize = 12;
                    fhUnits = fh.Units;
                    fh.Units = figParams.fig.units;
                    fhWH = fh.Position(3:4);
                    fh.Units = fhUnits;
                    txtStr = sprintf('You will need to scale the image X: %.2f  Y: %.2f to get appropriate dimensions.\nFor ref this fig should be %.1f x %.1f %s and this text font size should be %d\nPixel / point ratio is not 1.', ...
                        pxPtRatio(1),pxPtRatio(2), fhWH(1), fhWH(2), figParams.fig.units, mismatchedFontSize);
                    
                    % For vector formats, save on upper edge
                    % Because we assume these will be modified with editor
                    % like adobe illustrator.
                    if any(strcmpi(saveFormat,{'svg','pdf'}))
                        pos = [0 0.91 1 0.2];
                    else
                        % Otherwise, display at the top.
                        pos = [0 0.8 1 0.2];
                    end
                    
                    mismatchScale_h = annotation(fh,'textbox','LineStyle','none','String',txtStr, 'HorizontalAlignment','center',...
                        'Tag',scaleMsgTag, 'Color', 'r', ...
                        'FontSize', mismatchedFontSize, 'FontName', figParams.fig.fontName, ...
                        'Position', pos, 'Visible', 'on');
                end
        end
        
        function CloseFigIfOpen(figNum)
            fh = findobj('number', figNum);
            if ~isempty(fh)
                fh.Units = 'pixels';
                orgPos = fh.Position;
                close(fh)
                fh = figure(figNum);
                fh.Units = 'pixels';
                fh.Position = orgPos;
            end
            
            
        end
        
    end
    
    
    
    
end
    
