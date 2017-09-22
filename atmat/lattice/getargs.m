function varargout = getargs(args,default_values)
%GETARGS Process positional arguments from the input arguments
%
%ARGOUT=GETARGS(ARGIN,DEFAULT_VALUES) processes input arguments (typically
% from VARARGIN), by replacing DEFAULT_VALUES by valid ARGS items
% (items different from [], the empty numeric array)
%
%ARGOUT: all the arguments in a cell array, possibly longer than DEFAULT_VALUES:
%	length(ARGOUT)=max(length(ARGIN), length(DEFAULT_VALUES))
%
%[ARG1,ARG2,...]=GETARGS(ARGIN,DEFAULT_VALUES) returns arguments in
%	as many separate variables as element in DEFAULT_VALUES
%
%Example:
%
%function testfunc(varargin)
%
%[optflag,args]=getflag(varargin,'option');     % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
%[dname,dvalue]=getargs(args,{'abcd',297});     % Extract positional arguments
%
%See also GETFLAG, GETOPTION

na=length(default_values);
valid=~cellfun(@(arg) isempty(arg)&&isnumeric(arg),args);
default_values(valid)=args(valid);
if nargout==na
    varargout=default_values(1:na);
else
    varargout{1}=default_values;
end
end
