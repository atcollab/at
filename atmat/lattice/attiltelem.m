function elem = attiltelem(elem,rots)
%ATSHIFTELEM set new rotation parameters
%NEWELEM=ATROTELEM(OLDELEM,ROTS)
%
% ROTS - rotation angle in RADIANS
%   POSITIVE ROTS corresponds to a CORKSCREW (right) 
%   rotation of the ELEMENT looking in the direction of the beam.
%   (or CORKSCREW, aligned with s-axis) rotation of the ELEMENT
%   The rotation matrixes are stored in fields R1 and R2
%   R1 = [  cos(PSI) sin(PSI); -sin(PSI) cos(PSI) ]
%   R2 = R1'
%See also: atshiftelem, atmodelem

C=cos(rots);
S=sin(rots);
RM = diag([C C C C 1 1]);
RM(1,3) = S;
RM(2,4) = S;
RM(3,1) = -S;
RM(4,2) = -S;
elem.R1=RM;
elem.R2=RM';
end

