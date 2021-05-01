function [rsrc,varargout] = decodeatargs(default_values,args)
%DECODEATARGS separates arguments and resources
%
%  [RSRC,ARGS]=decodeatargs(DEFARGS,ARGLIST)
%
%  INPUTS
%    1. DEFARGS - Default values of mandatory argument
%    2. ARGS    - Arguments
%
%  OUPUTS
%    1. rsrc      - Optional arguments
%    2. varargout - Mandatory arguments
%
%  See also getoption, getflag

na=length(default_values);
mandatory = @(arg) ~ischar(arg) || isempty(arg) || endsWith(arg, 'Pass');
[varargout{1:na},rsrc]=getargs(args,default_values{:},'check',mandatory);
end
