function NEWPARAMGROUP = saveparamgroup(LATTICE,PARAMGROUP)
%SAVEPARAMGROUP saves the values of multiple physical 
% parameters of the lattice in the special SavedValue field of 
% AT parameter group structure. The result can be late used 
% with RESTOREPARAMGROUP
%
% PARAMGROUP = saveparamgroup(LATTICE,PARAMGROUP)
% 
% See also: ATPARAMGROUP RESTORPARAMGROUP SETPARAMGROUP

NEWPARAMGROUP=PARAMGROUP;

for i=1:length(PARAMGROUP)
   NEWPARAMGROUP(i).SavedValue = ...
       getfield(LATTICE{PARAMGROUP(i).ElemIndex},...
       PARAMGROUP(i).FieldName,PARAMGROUP(i).FieldIndex);
end
