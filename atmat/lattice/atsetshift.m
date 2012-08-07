function ring=atsetshift(varargin)
%ATSETTILT sets the misalignment vectors
%
% RING=ATSETSHIFT(RING, ELEMINDEX, DX, DY) sets the entrance and exit misalignment vectors
%  of one element or a group of elements in the globally defined lattice THERING.
%
% ELEMINDEX contains indexes of elements to be misaligned
% DX, DY are displacements of the ELEMENT
%  so the coordinate transformation on the particle at entrance is
%	X  ->  X-DX
%   Y  ->  Y-DY
%
% ATSETSHIFT(ELEMINDEX, DX, DY) Uses the global variable THERING
%
% See also: ATSETTILT

global THERING
if ~iscell(varargin{1})
    THERING=atsetshift(THERING,varargin{:});
else
    [ring,idx,dx,dy]=deal(varargin{:});
    
    if length(dx) == 1
        dx=dx*ones(size(idx));
    elseif length(dx) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(dx));
    end
    if length(dy) == 1
        dy=dy*ones(size(idx));
    elseif length(dy) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(dy))
    end
    
    for i = 1:length(idx)
        ring{idx(i)}=atshiftelem(ring{idx(i)},dx(i),dy(i));
    end
end
