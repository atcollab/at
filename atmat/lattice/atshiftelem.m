function elem = atshiftelem(elem,dx,dz)
%ATSHIFTELEM set new displacement parameters
%NEWELEM=ATSHIFTELEM(OLDELEM,DX,DZ)
%
%DX:	horizontal element displacement
%DZ:	Vertical element displacement
%
%See also: attiltelem, atmodelem

disp=[dx;0;dz;0;0;0];
elem.T1=-disp;
elem.T2=disp;
end
