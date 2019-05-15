function elem = attiltelem(elem,PSI,varargin)
%ATTILTELEM sets new rotation parameters
%  NEWELEM = ATTILTELEM(OLDELEM,PSI)
%
%  PSI - rotation angle in radians
%   POSITIVE ROTS corresponds to a CORKSCREW (right) 
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The rotation matrixes are stored in fields R1 and R2
%   R1 = [cos(PSI) sin(PSI); -sin(PSI) cos(PSI)]
%   R2 = R1'
%
%   NEWELEM=ATTILTELEM(OLDELEM,PSI,'RelativeTilt')
%   the rotation is added to the previous misalignment
%
%  See also ATSETTILT, ATSHIFTELEM, ATMODELEM

[RelativeTilt,~]=getflag(varargin,'RelativeTilt');
C       = cos(PSI);
S       = sin(PSI);
RM      = diag([C C C C 1 1]);
RM(1,3) = S;
RM(2,4) = S;
RM(3,1) = -S;
RM(4,2) = -S;
if RelativeTilt && isfield(elem,'R1') && isfield(elem,'R2')
    newR1=RM*elem.R1;
    newR2=RM'*elem.R2;
    elem.R1=newR1;
    elem.R2=newR2;
else
    elem.R1 = RM;
    elem.R2 = RM';
end
end

