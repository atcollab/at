function [opts,optval] = getoption(opts,optname,optval)
%GETOPTION Extract one option from an option list ['Name',value,...]
%
%[NEWOPTLIST,OPTVAL]=GETOPTION(OPTLIST,OPTNAME,OPTDEFAULT)
%
%OPTLIST:   Option list
%OPTNAME:   Name of the desired option
%OPTDEFAULT:Default value for the option (defaul: [])
%
ok=strcmpi(optname,opts(1:2:end));
if any(ok)
    optval=opts{2*find(ok,1,'last')};
    opts(reshape([ok;ok],1,[]))=[];
end
end

