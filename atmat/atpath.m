function atpath(varargin)
%ATPATH adds the AT directories to the MATLAB search path
    
ATROOT = fileparts(mfilename('fullpath'));
ATINTEGRATORS = fullfile(fileparts(ATROOT),'atintegrators');
disp(ATROOT)
disp(ATINTEGRATORS)
addpath(genpath(ATROOT));
addpath(genpath(ATINTEGRATORS));

