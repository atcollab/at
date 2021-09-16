function [linargs,opts] = linoptions(opts, dp)
%LINOPTIONS  (private) extract arguments for atlinopt
%
% convert an argument list like
% (0.01,'opt1', val1) to
% ('dp',0.01,'opt1',val1)
%or
% (0, 'opt1', val1) to
% ('opt1',val1)

[linargs,opts]=getoption(opts,{'dp','dct','twiss_in'});
if nargin >=2 && dp ~= 0
    linargs=[{'dp',dp} linargs];
end