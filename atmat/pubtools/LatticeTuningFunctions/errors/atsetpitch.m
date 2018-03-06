function ring=atsetpitch(varargin)
%ATSETPITCH sets the entrance and exit rotation matrices
% of an element or a group of elements in THERING
%
% RING=ATSETPITCH(RING,ELEMINDEX, PSI)
% ELEMINDEX contains indexes of elements to be rotated
% PHI - angle(s) of rotation in RADIANS
%   POSITIVE PHI corresponds to a CORKSCREW (down)
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The misalgnment matrixes are stored in fields T1 and T2 adding to the
%   existing T1, T2
%
% ATSETPITCH(ELEMINDEX, PHI) Uses the global variable THERING
%
% See also ATSETYAW

global THERING
if ~iscell(varargin{1})
    THERING=atsetpitch(THERING,varargin{:});
else
    [ring,idx,rot]=deal(varargin{:});
    
    if length(rot) == 1
        rot=rot*ones(size(idx));
    elseif length(rot) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rot))
    end
    
    for i = 1:length(idx)
        ring{idx(i)}=pitchelem(ring{idx(i)},rot(i));
    end
end

end

function elem = pitchelem(elem,PHI)
%ATTILTELEM sets new rotation parameters
%  NEWELEM = pitchelem(OLDELEM,PSI)
%
%  PHI - rotation angle in radians
%   POSITIVE ROTS corresponds to a CORKSCREW (down) 
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The transformation vectors are stored in fields T1 and T2
%   T1 = T1(existing) + [0 0 -L/2*sin(PSI) PSI 0 0]'
%   T2 = T2(existing) + [0 0 -L/2*sin(PSI) - PSI 0 0]'
% 
%  See also ATSHIFTELEM, ATMODELEM

if isfield(elem,'Length')
    C   = sin(PHI);
    Lh  = elem.Length/2;
    
    T10=zeros(6,1);
    T20=T10;
    if isfield(elem,'T1')
        T10 = elem.T1;
        T20 = elem.T2;
    end
    
    elem.T1 = T10 + [0; 0; -Lh*C; PHI; 0; 0];
    elem.T2 = T20 + [0; 0; -Lh*C; -PHI; 0; 0];
end

end
