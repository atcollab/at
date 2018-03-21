function orbit  = findorbit(RING,D, varargin)
%FINDORBIT Alias to the orbit search functions: (findorbit4, findsyncorbit,
%          findorbit4 ...) 
%  depending on the number of parameters, their types and argument options
%  It can also return additional parameters.
%  
%  NOTES
%  1. Temporary version of 7/18/00 - defaults to FINDORBIT4
%
% See also FINDORBIT4, FINDSYNCORBIT, FINDORBIT6.

orbit = findorbit4(RING, D, varargin{:});

return