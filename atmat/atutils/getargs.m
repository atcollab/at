function varargout = getargs(args,varargin)
%GETARGS Process positional arguments from the input arguments
%
%Processes input arguments (typically from VARARGIN), by replacing default values
%by valid ARGIN items (items different from [], the empty numeric array)
%
%[V1,V2,...,REMARGS]=GETARGS(ARGIN,DEF1,DEF2,...) returns as many variables
% as default values plus a cell array of remaining arguments.
%
%Example:
%
%function testfunc(varargin)
%
%[optflag,args]=getflag(varargin,'option');     % Extract an optional flag
%[range,args]=getoption(args,'Range', 1:10);	% Extract a keyword argument
%[dname,dvalue]=getargs(args,'abcd',297);     % Extract positional arguments
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
