function indexstruct = atindex(lattice,varargin)
%ATINDEX extracts the information about element families and
% indexes from AT lattice
% ATI = ATINDEX(LATTICE)
%  returns a srtucture with fields named after element family
%  containing an array of their AT indexes;
%
%   ATI.QF = [...]
%   ATI.QD = [...];
%   ...
% See also FINDCELLS

indexstruct=struct();
for i = 1:length(lattice);
    if isfield(lattice{i},'FamName') && ~isempty(lattice{i}.FamName)
        famname=lattice{i}.FamName;
        try
            a.(famname)=0; %#ok<STRNU>
        catch %#ok<CTCH>
            famname=['x' famname]; %#ok<AGROW>
        end
    else
        famname='unnamed';
    end
    if isfield(indexstruct,famname)
        indexstruct.(famname)(end+1)=i;
    else
        indexstruct.(famname)(1)=i;
    end
end
