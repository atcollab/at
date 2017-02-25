function atsetglobval(Energy)
%ATSETGLOBVAL creates the global variable GLOBVAL and adds Energy
%as a field
%
%  NOTES
%    1. This use of global variables is deprecated, but for backwards
%       compatibility, this is sometimes necessary.
%
%  See also ATDISPLAY

global GLOBVAL

GLOBVAL.E0 = Energy;