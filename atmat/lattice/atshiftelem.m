function elem = atshiftelem(elem,dx,dy,varargin)
%ATSHIFTELEM set new displacement parameters
%NEWELEM=ATSHIFTELEM(OLDELEM,DX,DZ)
%
%DX:	horizontal element displacement
%DY:	Vertical element displacement
%
%
%NEWELEM=ATSHIFTELEM(OLDELEM,DX,DY,'RelativeShift')
%the shift is added to the previous misalignment of the element
%
%See also: atsetshift, attiltelem, atmodelem

RelativeShift=getflag(varargin,'RelativeShift');
disp=[dx;0;dy;0;0;0];
if RelativeShift && isfield(elem,'T1') && isfield(elem,'T2')
    elem.T1=elem.T1-disp;
    elem.T2=elem.T2+disp;
else
    elem.T1=-disp;
    elem.T2=disp;
end
end