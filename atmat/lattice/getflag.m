function [flag,opts] = getflag(opts,optname)
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
%Dee also GETOPTION

ok=strcmpi(optname,opts);
flag=any(ok);
opts=opts(~ok);
end

