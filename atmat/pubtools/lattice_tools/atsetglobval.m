function  atsetglobval(Energy)
%atsetglobval(Energy) creates the global variable GLOBVAL and adds Energy
%as a field
%This use of global variables is deprecated, but for backwards
%compatibility, this is sometimes necessary.

global GLOBVAL

GLOBVAL.E0=Energy;