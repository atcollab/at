function [flag,opts] = getflag(opts,optname)
%GETFLAG Check the presence of an option in an argument list
%
%OPTION=GETFLAG(ARGS,OPTNAME)
%
%ARGS:      Argument list (cell array)
%OPTNAME:	Name of the desired option (string)
%
%[OPTION,NEWARGS]=GETFLAG(ARGS,OPTNAME)
%           Returns the argument list after removing the processed flag

ok=strcmpi(optname,opts);
flag=any(ok);
opts=opts(~ok);
end

