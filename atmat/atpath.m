function atpath(varargin)
%ATPATH adds the AT directories to the MATLAB search path
    
ATROOT = fileparts(mfilename('fullpath'));
ATPATH = genpath(ATROOT);
disp(ATROOT);
addpath(ATPATH);
