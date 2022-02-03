function [linargs,opts] = linoptions(opts)
%LINOPTIONS  (private) extract arguments for atlinopt
%
% Separate the options for atlinopt

[linargs,opts]=getoption(opts,{'dp','dct','df','twiss_in','orbit',...
    'XYStep','DPStep','OrbConvergence','OrbMaxIter','guess','method'});
