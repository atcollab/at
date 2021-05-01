function varargout = parseargs(default_values,args)
%PARSEARGS Check and expands optional argument lists
%ARGOUT=PARSEARGS(DEFAULT_VALUES,ARGIN)
%[ARG1,ARG2,...]=PARSEARGS(DEFAULT_VALUES,ARGIN)
%
% obsolete: see GETARGS

[varargout{1:nargout}]=getargs(args,default_values{:});
end
