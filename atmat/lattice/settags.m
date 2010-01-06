function LATTICE = settags(LATTICE,INDEX,tag,varargin)
%SETTAGS sets the 'Tag' field in AT lattice elements
% LATTICE = SETTAGS(LATTICE, INDEX, TAG) 
%   INDEX can be integer AT index or a string famly name
%   TAG is a string tag or a cell array of strings
% LATTICE = SETTAGS(LATTICE, INDEX, TAG, 'append') 
%   appends to existing tags

if ischar(INDEX)
    INDEX = findcells(LATTICE,'FamName',INDEX);
    if isempty(INDEX)
        error(['Family ''',fname,''' is not found']);  
    end
end

tag = reshape(cellstr(tag),1,[]);

if nargin <=3
    for i = INDEX
        LATTICE{i}.Tag = tag;
    end
elseif strcmpi(varargin{1},'append')
    for i = INDEX
        if isfield(LATTICE{i},'Tag')
            LATTICE{i}.Tag = reshape(cellstr(LATTICE{i}.Tag),1,[]);
            LATTICE{i}.Tag = union(LATTICE{i}.Tag, tag);
        else
            LATTICE{i}.Tag = tag;
        end
        
    end
end
