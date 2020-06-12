function varargout = getargs(args,default,varargin)
%GETARGS Process positional arguments from the input arguments
%
%VALUES=GETARGS(ARGIN,DEFAULT_VALUES) processes input arguments (typically
% from VARARGIN), by replacing DEFAULT_VALUES by valid ARGIN items
% (items different from [], the empty numeric array)
%
%VALUES: all the arguments in a cell array, possibly longer than DEFAULT_VALUES:
%	length(ARGOUT)=max(length(ARGIN), length(DEFAULT_VALUES))
%
%[V1,V2,...,REMAIN]=GETARGS(ARGIN,DEFAULT_VALUES) returns arguments in
%	as many separate variables as element in DEFAULT_VALUES and adds a cell
%   array of remaining arguments, if any
%
%[...]=GETARGS(ARGIN,DEF1,DEF2,...) takes default values in separate arguments
%
%Example:
%
%function testfunc(varargin)
%
%[optflag,args]=getflag(varargin,'option');     % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
%[dname,dvalue]=getargs(args,{'abcd',297});     % Extract positional arguments
%or
%[dname,dvalue]=getargs(args,'abcd',297);     % Extract positional arguments
%
%See also GETFLAG, GETOPTION

default_values=[default varargin];
na=length(default_values);
valid=~cellfun(@(arg) isempty(arg)&&isnumeric(arg),args);
default_values(valid)=args(valid);
if nargout>=na
    varargout=[default_values(1:na) {default_values(na+1:end)}];
else
    varargout{1}=default_values;
end
end
