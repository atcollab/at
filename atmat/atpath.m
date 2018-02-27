function atpath(varargin)
%ATPATH Adds the AT directories to the MATLAB search path
    
ATROOT = fileparts(mfilename('fullpath'));
ATINTEGRATORS = fullfile(fileparts(ATROOT),'atintegrators');
MACHINEDATA = fullfile(fileparts(ATROOT),'machine_data');
disp(ATROOT)
addpath(genpath(ATROOT));
disp(ATINTEGRATORS)
addpath(genpath(ATINTEGRATORS));
disp(MACHINEDATA)
addpath(genpath(MACHINEDATA));