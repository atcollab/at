function [optval,opts] = getoption(opts,optname,defval)
%GETOPTION Extract a keyword argument from the input arguments
%
%OPTVAL=GETOPTION(ARGS,OPTNAME)
%   Extract a keyword argument, in the form of a pair "name,value" from
%   input arguments (typically from VARARGIN).
%	Return the value of the desired option from the argument list and
%	send an exception if it is absent. Option names are case insensitive
% ARGS:     Argument list (cell array or structure)
% OPTNAME:	Name of the desired option
%
%
%OPTVAL=GETOPTION(ARGS,OPTNAME,OPTDEFAULT)
%	Return a default value if the option is absent
% OPTDEFAULT:   Default value for the option
%
%[OPTVAL,NEWARGS]=GETOPTION(...)
%  Return the remaining arguments after removing the processed ones
%
%Example:
%
%function testfunc(varargin)
%
%[optflag,args]=getflag(varargin,'option');     % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
%[width, height]=getargs(args,{210,297});       % Extract positional arguments
%
%See also GETFLAG, GETARGS

if iscell(opts)
    ok=[strcmpi(optname,opts(1:end-1)) false];  %option name cannot be the last argument
    if any(ok)
        okval=circshift(ok,[0,1]);
        defval=opts{find(okval,1,'last')};
        opts(ok|okval)=[];
    end
elseif isstruct(opts)
    if isfield(opts,optname)
        defval=opts.(optname);
        opts=rmfield(opts,optname);
    end
end
try
    optval=defval;
catch       % Dould be 'MATLAB:minrhs'
    error('getoption:missing','Option "%s" is missing',optname);
end
end
