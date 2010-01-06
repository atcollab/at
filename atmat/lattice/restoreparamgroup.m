function NEWLATTICE = restoreparamgroup(LATTICE,PARAMGROUP)
%RESTOREPARAMGROUP restores the values of multiple physical 
% parameters of the lattice. 
% NEWLATTICE = RESTOREPARAMGROUP(LATTICE,PARAMGROUP)
% 
% See also: ATPARAMGROUP RESTORPARAMGROUP SAVEPARAMGROUP
NEWLATTICE=LATTICE;

for i=1:length(PARAMGROUP)
   NEWLATTICE{PARAMGROUP(i).ElemIndex}=...
       setfield(NEWLATTICE{PARAMGROUP(i).ElemIndex},...
       PARAMGROUP(i).FieldName,PARAMGROUP(i).FieldIndex,...
       PARAMGROUP(i).SavedValue);
end
