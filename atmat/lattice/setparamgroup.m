function NEWLATTICE = setparamgroup(LATTICE,PARAMGROUP,PVALUE,varargin)
%SETPARAMGROUP modifies a group of parameters
% NEWLATTICE = setparamgroup(LATTICE,PARAMGROUP,PVALUE)
% 
% See also: ATPARAMGROUP RESTORPARAMGROUP

NEWLATTICE=LATTICE;


if nargin == 3
    for i=1:length(PARAMGROUP)
        NEWLATTICE{PARAMGROUP(i).ElemIndex}=...
            setfield(NEWLATTICE{PARAMGROUP(i).ElemIndex},...
            PARAMGROUP(i).FieldName,PARAMGROUP(i).FieldIndex,...
            feval(PARAMGROUP(i).Function,PVALUE,PARAMGROUP(i).Args{:}));
    end
else
    for i=1:length(PARAMGROUP)
        NEWLATTICE{PARAMGROUP(i).ElemIndex}=...
            setfield(NEWLATTICE{PARAMGROUP(i).ElemIndex},...
            PARAMGROUP(i).FieldName,PARAMGROUP(i).FieldIndex,...
            feval(PARAMGROUP(i).Function,PVALUE,varargin{:}));
    end
end
