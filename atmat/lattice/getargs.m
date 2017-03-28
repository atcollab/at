function varargout = getargs(args,default_values)
%GETARGS Check and expands optional argument lists
%ARGOUT=GETARGS(ARGIN,DEFAULT_VALUES)
%[ARG1,ARG2,...]=GETARGS(ARGIN,DEFAULT_VALUES)

na=min(length(default_values),length(args));
valid=~cellfun(@(arg) isempty(arg)&&isnumeric(arg),args(1:na));
default_values(valid)=args(valid);
if nargout==length(default_values)
    varargout=default_values;
else
    varargout{1}=default_values;
end
end
