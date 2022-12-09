function [rsrc,varargout] = decodeatargs(default_values,args)
%DECODEATARGS separates arguments and resources
%
%  [RSRC,ARGS]=decodeatargs(DEFARGS,ARGLIST)
%
%  INPUTS
%    1. DEFARGS - Default values for mandatory argument
%    2. ARGLIST - Arguments
%
%  OUPUTS
%    1. RSRC    - Optional arguments (all remaining arguments)
%    2. ARGS    - Mandatory arguments
%
%  See also getoption, getflag

na=length(default_values);
mandatory = @(arg) ~(ischar(arg) || isstring(arg)) || endsWith(arg, 'Pass');
[varargout{1:na},rsrc]=getargs(args,default_values{:},'check',mandatory);
end
