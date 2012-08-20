function z=rbend(fname,L,A,A1,A2,K,method)
%BEND('FAMILYNAME',  Length[m], BendingAngle[rad], EntranceAngle[rad],
%	ExitAngle[rad], K, 'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName        	family name
%		Length         	length of the arc for an on-energy particle [m]
%		BendingAngle		total bending angle [rad]
%		EntranceAngle		[rad] (A/2 - for rectangular bends)
%		ExitAngle			[rad] (A/2 - for rectangular bends)
%		ByError				error in the dipole field relative to the design value 
%		K						quadrupole K-value for combined funtion bends
%		PassMethod        name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

ElemData=atrbend(fname,L,A,K,method);
Elem.EntranceAngle=A1;  %for backwards compatibility
Elem.ExitAngle=A2;

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

