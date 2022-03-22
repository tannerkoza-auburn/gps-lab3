%% startup.m

% DESCRIPTION: startup.m is a file that is ran when MATLAB is started with
% this working directory. Any desired commands can be ran during startup.
% If switching to a new directory that contains a startup.m file in an open 
% instance of MATLAB, type startup in the Command Window to call desired 
% startup commands.

%% Functionality

% Add All Subfolders of Current Folder
dir = fileparts(which(mfilename));
addpath(genpath(dir))