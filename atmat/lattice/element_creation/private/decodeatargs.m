function [rsrc,varargout] = decodeatargs(default_values,args)
%DECODEATARGS separates arguments and resources
%
%  [RSRC,ARGS]=decodeatargs(DEFARGS,ARGLIST)
%
%  INPUTS
%    1. DEFARGS - Values per default if not specify for mandatory
%                 argument
%    2. ARGS    - Other arugments
%
%  OUPUTS
%    1. rsrc      - Mandatory arguments
%    2. varargout - Optional arguments
%
%  NOTES
%    1. DEFARGS must have length >= 2
%
%  See also getoption, getflag

% function to get the first string which is after the PassMethod
nopass    = @(arg) ischar(arg) && ~isempty(arg) && isempty(regexp(arg,'.*Pass$','once'));
% get index of first optional char Argument
chararg   = find(cellfun(nopass,[args {'x'}]),1); 
% get all optional argument (ressources)
rsrc      = args(chararg:end);
% output mandatory argument
varargout = parseargs(default_values,args(1:chararg-1));

end
