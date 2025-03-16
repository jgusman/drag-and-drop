
%% Necessary inputs

%% Necessary outputs

% Neural data can live as a class property

classdef Decoder < matlab.mixin.Copyable
    properties
        % Required properties
        name % Name of the decoder
        mode % Mode - discrete, kinematic, idle
%         params % Decoder params
        
        % Suggested properties
%         data  
%         labels 
%         coefficients
%         info
%         performance
    end
    methods
        
        function Train(obj)
            % Train/Calibrate the decoder
        
        end
        
        function Test(obj)
            % Test your data
        
        end
        
        function Summarize(obj)
            % Plot/summarize/save results
        
        end
    end
    
    %% Helper functions 
    methods
        % -----------------------------------------------------------------
        %   Save
        % -----------------------------------------------------------------
        function sobj = saveobj(obj)
        % Save property fields, but do not save the object.
        % TO DO: have a function that can load this struct to recreate the
        % decoder object
             fn = fieldnames(obj);
             for ii = 1:length(fn)
                 sobj.(fn) = obj.(fn);
             end

             if ~obj.params.onSave.keepDataAndLabels
                 sobj.data = [];
                 sobj.labels = [];
             end    
        end
        
        
        
    end
end