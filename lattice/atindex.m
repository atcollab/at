function indexstruct = atindex(LATTICE,varargin);
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
indexstruct = struct([]);
fieldnames = {};
indexes = {};
for i = 1:length(LATTICE);
    
    if isfield(LATTICE{i},'FamName') & ~isempty(LATTICE{i}.FamName);
        m = find(strcmp(fieldnames,LATTICE{i}.FamName));
        if isempty(m)
            fieldnames{end+1} = LATTICE{i}.FamName;
            indexes{end+1}{1} = i;
        else
            indexes{m}{end+1}=i; 
        end
    end
end

for m = 1:length(indexes)
    indexes{m}=cell2mat(indexes{m});
end
ARGS = reshape([fieldnames;indexes],1,[]);
indexstruct = struct(ARGS{:});