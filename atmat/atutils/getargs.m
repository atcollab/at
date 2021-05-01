function varargout = getargs(args,varargin)
%GETARGS Process positional arguments from the input arguments
%
%Process input arguments (typically from VARARGIN), replacing missing
%arguments by default values. The default values is also used when an
%argument is "[]" (empty numeric array).
%
%[V1,V2,...,REMARGS]=GETARGS(ARGIN,DEF1,DEF2,...)
%   Return as many variables as default values, plus a cell array of
%   remaining arguments.
%
%[...]=GETARGS(ARGIN,...,'check',checkfun)
%   Check each input arguments using checkfun(arg) and stops processing
%   at the first failing argument. Remaining arguments are available in
%   REMARGS.
%
%Example 1:
%
%[optflag,args]=getflag(varargin,'option');   % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);  % Extract a keyword argument
%[dname,dvalue]=getargs(args,'abcd',297);     % Extract positional arguments
%
%Example 2:
%
%global THERING
%[ring,args]=getargs(varargin,THERING,'check',@iscell)
%   if the 1st argument is a cell array, it will be used as "ring",
%   otherwise, THERING will be used. In any case, the remaining
%   arguments are availble in "args".
%
%See also GETFLAG, GETOPTION

[check,default_values]=getoption(varargin,'check',@(arg) true);
na=length(default_values);
% Look for default args
takedef = cellfun(@(arg) isempty(arg) && isnumeric(arg),args);
% Look for valid arguments and stop at 1st non-valid
checked =  takedef | cellfun(@(arg) check(arg),args);
checked(find(~checked,1):end)=false;
% Look for valid, non-default arguments
valid = checked & ~takedef;
default_values(valid)=args(valid);
% Append non-valid arguments
default_values=[default_values args(~checked)];
varargout=[default_values(1:na) {default_values(na+1:end)}];
end
