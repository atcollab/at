function varargout = parseargs(default_values,args)
%PARSEARGS Check and expands optional argument lists
%ARGOUT=PARSEARGS(DEFAULT_VALUES,ARGIN)
%[ARG1,ARG2,...]=PARSEARGS(DEFAULT_VALUES,ARSIN)

na=min(length(default_values),length(args));
ok=~cellfun(@isempty,args(1:na));
default_values(ok)=args(ok);
if nargout==length(default_values)
    varargout=default_values;
else
    varargout{1}=default_values;
end
