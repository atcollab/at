function [dst,src] = mvfield(dst,src,fields)
%MVFIELD Move fields from one structure to another
%
%[NEWDST,NEWSRC]=MVFIELD(DST,SRC,FIELDNAMES)
%   Moves the selected fields from SRC to DST
%
%DST:           Destination structure
%SRC:           Source structure
%FIELDNAMES:    Field names to be moved (Default: all fields from SRC)
%
%NEWDST:        DST structure with fields added
%NEWSRC:        SRC structure with fields removed

if nargin >= 3      % Select valis fields
    fnames=fields(isfield(src,fields));
else                % Take all fields from src
    fnames=fieldnames(src);
end
cellfun(@movefields,fnames);

    function movefields(fld)
        dst.(fld)=src.(fld);
        src=rmfield(src,fld);
    end
end

