function ring=atsetyaw(varargin)
%ATSETYAW sets the entrance and exit rotation matrices
% of an element or a group of elements in THERING
%
% RING=ATSETYAW(RING,ELEMINDEX, PSI)
% ELEMINDEX contains indexes of elements to be rotated
% THETA - angle(s) of rotation in RADIANS
%   POSITIVE THETA corresponds to a CORKSCREW (left)
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The misalgnment matrixes are stored in fields T1 and T2 adding to the
%   existing T1, T2
%
% ATSETYAW(ELEMINDEX, THETA) Uses the global variable THERING
%
% See also ATSETPITCH

global THERING
if ~iscell(varargin{1})
    THERING=atsetyaw(THERING,varargin{:});
else
    [ring,idx,rot]=deal(varargin{:});
    
    if length(rot) == 1
        rot=rot*ones(size(idx));
    elseif length(rot) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rot))
    end
    
    for i = 1:length(idx)
        ring{idx(i)}=yawelem(ring{idx(i)},rot(i));
    end
end

end

function elem = yawelem(elem,THETA)
%ATTILTELEM sets new rotation parameters
%  NEWELEM = pitchelem(OLDELEM,THETA)
%
%  THETA - rotation angle in radians
%   POSITIVE THETA corresponds to a CORKSCREW (left) 
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The transformation vectors are stored in fields T1 and T2
%   T1 = T1(existing) + [L/2*sin(THETA) THETA 0 0 0 0]'
%   T2 = T2(existing) + [L/2*sin(THETA) - THETA 0 0 0 0]'
% 
%  See also ATSHIFTELEM, ATMODELEM

if isfield(elem,'Length')
    C   = sin(THETA);
    Lh  = elem.Length/2;
    
    T10=zeros(6,1);
    T20=T10;
    if isfield(elem,'T1')
        T10 = elem.T1;
        T20 = elem.T2;
    end
    
    elem.T1 = T10 + [-Lh*C; THETA; 0; 0; 0; 0];
    elem.T2 = T20 + [-Lh*C; -THETA; 0; 0; 0; 0];
end

end
