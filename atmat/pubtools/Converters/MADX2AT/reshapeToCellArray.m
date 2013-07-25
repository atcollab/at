function CEL_AR=reshapeToCellArray(CEL_CEL)
% if CEL_CEL is a cell array of structures and cell arrays it converts it a
% cell array of structures.

CEL_AR={};

for i=1:length(CEL_CEL)
    if ~isstruct(CEL_CEL{i})
      CEL_AR=[CEL_AR; reshapeToCellArray(CEL_CEL{i})]; %#ok<AGROW>
    else
      CEL_AR=[CEL_AR; CEL_CEL{i}];
    end
end


