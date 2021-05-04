function [flag,opts] = getflag(opts,key)
%GETFLAG Check the presence of a flag in an argument list
%
%OPTION=GETFLAG(ARGS,OPTNAME)
%   Return a logical value indicating the presence of the flag name in the
%   argument list. Flag names are case insensitive.
%
%ARGS:      Argument list (cell array)
%OPTNAME:	Name of the desired option (string)
%
%[OPTION,NEWARGS]=GETFLAG(ARGS,OPTNAME)
%           Returns the argument list after removing the processed flag
%
%Example:
%
%function testfunc(varargin)
%
%[optflag,args]=getflag(varargin,'option');     % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
%[width, height]=getargs(args, 210, 297);       % Extract positional arguments
%
%Dee also GETOPTION, GETARGS

ok = cellfun(@(v) (ischar(v) || isstring(v)) && strcmpi(v, key),opts);
flag = any(ok);
opts=opts(~ok);
end
