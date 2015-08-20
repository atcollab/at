function [optval,opts] = getoption(opts,optname,optval)
%GETOPTION Extract one option from an option list ['Name',value,...]
%
%OPTVAL=GETOPTION(OPTIONS,OPTNAME,OPTDEFAULT)
%
%OPTLIST:   Option list (cell array or structure)
%OPTNAME:   Name of the desired option
%OPTDEFAULT:Default value for the option (default: [])
%
%[OPTVAL,NEWOPTIONS]=GETOPTION(...)
%           Returns new options after removing the processed ones
%
%See also GETFLAG

if iscell(opts)
    ok=strcmpi(optname,opts(1:2:end));
    if any(ok)
        optval=opts{2*find(ok,1,'last')};
        opts(reshape([ok;ok],1,[]))=[];
    end
elseif isstruct(opts) && isfield(opts,optname)
    optval=opts.(optname);
    opts=rmfield(opts,optname);
end
end

