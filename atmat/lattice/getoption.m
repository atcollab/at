function [optval,opts] = getoption(opts,optname,defval)
%GETOPTION Extract one option from an option list ['Name',value,...]
%
%OPTVAL=GETOPTION(ARGS,OPTNAME)
%	Return the value of the desired option from the argument list and
%	send an exception if it is absent. Option names are case insensitive
% OPTLIST:      Option list (cell array or structure)
% OPTNAME:      Name of the desired option
%
%
%OPTVAL=GETOPTION(ARGS,OPTNAME,OPTDEFAULT)
%	Return a default value if the option is absent
% OPTDEFAULT:   Default value for the option
%
%[OPTVAL,NEWARGS]=GETOPTION(...)
%  Returns remaining options after removing the processed ones
%
%See also GETFLAG

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
catch
    error('getoption:missing','Option "%s" is missing',optname);
end
end
