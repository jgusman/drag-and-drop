% Main script for generating paper figures
%
% Code released with manuscript: Gusman, Hosman, et al., "Multi-gesture drag-and-drop
% decoding in a 2D iBCI control task", Journal of Neural Engineering (2025).
%
% Copyright Jacob Gusman, 2025. Brown University
% jacob_gusman@brown.edu
% -------------------------------------------------------------------------


%% Add paths and specify outputs and data directories
activeFilePath = matlab.desktop.editor.getActiveFilename; % gets absolute path to the current file (DragAndDropMaster.m)
scriptDir = fileparts(activeFilePath);                    % gets the scripts folder path
repoDir = fileparts(scriptDir);                           % gets the repo path
addpath(genpath(fullfile(repoDir)));                      % adds folder and subfolders to matlab path

saveFiguresFolder = fullfile(repoDir,'outputs');          % set folder in with saved generate figures (defaults to 'outputs' folder)
dataFolder = '';                                          % set folder where the data from Dryad was saved locally




%% %%%%%%%%%%%%%%%%%%%% Load Gesture Hero task data %%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(dataFolder,'sesData_GH_T11.mat'))
load(fullfile(dataFolder,'sesData_GH_T5.mat'))

%% KW sweep on Gesture Hero (Figs 2A, 3A-B, S2A, S3A, S4A, S5)
GH_SelectivityOverTime

%% LDA over time (Figs 2B-C, S2B-C, S3B-C, S4B-C)
GH_LDAOverTime

%% Neural components grid search (Figs 3C-D,S6, S7, S8)
GH_NeuralComponentsSearch

%% Latch Decoder example trial and comparison with Gesture Decoder (Figs 4B-D)
GH_LatchDecoder

%% Feature selection analysis (Fig 5)
GH_FeatureSelection



%% %%%%%%%%%%%%%%%%%%%% Load Drag and Drop task data %%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(dataFolder,'sesData_DD_T11.mat'))

%% Time to target plots (Figs 6C-D)
DD_TimeToTarget

%% Gesture Decoder (offline) vs. Latch Decoder (Fig 7)
DD_GestureVsLatch

