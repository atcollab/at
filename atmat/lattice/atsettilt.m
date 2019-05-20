function ring=atsettilt(varargin)
%ATSETTILT sets the entrance and exit rotation matrices
% of an element or a group of elements in THERING
%
% RING=ATSETTILT(RING,ELEMINDEX, PSI)
% ELEMINDEX contains indexes of elements to be rotated
% PSI - angle(s) of rotation in RADIANS
%   POSITIVE PSI corresponds to a CORKSCREW (right)
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The misalgnment matrixes are stored in fields R1 and R2
%   R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
%   R2 = R1'
%
% RING=ATSETTILT(RING,ELEMINDEX,PSI,'RelativeTilt')
% the rotation is added to the previous misalignment
%
% ATSETTILT(ELEMINDEX, PSI) Uses the global variable THERING
%
% See also ATSETSHIFT

global THERING
if ~iscell(varargin{1})
    THERING=atsettilt(THERING,varargin{:});
else
    [ring,idx,rot]=deal(varargin{1:3});
    
    if length(rot) == 1
        rot=rot*ones(size(idx));
    elseif length(rot) ~= length(idx)
        error('AT:length','Vector lengths are incompatible: %i/%i.',length(idx),length(rot))
    end
    for i = 1:length(idx)
        ring{idx(i)}=attiltelem(ring{idx(i)},rot(i),varargin{4:end});
    end
end
