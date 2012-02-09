function Elem=atsbend(fname,L,A,A1,A2,K,method)
%BEND('FAMILYNAME',  Length[m], BendingAngle[rad], EntranceAngle[rad],
%	ExitAngle[rad], K, 'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName        	family name
%		Length         	length of the arc for an on-energy particle [m]
%		BendingAngle		total bending angle [rad]
%		EntranceAngle		[rad] (0 - for sector bends)
%		ExitAngle			[rad] (0 - for sector bends)
%		ByError				error in the dipole field relative to the design value 
%		K						quadrupole K-value for combined funtion bends
%		PassMethod        name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family


Elem.FamName = fname;  % add check for identical family names
Elem.Length			= L;
Elem.MaxOrder			= 3;
Elem.NumIntSteps 	= 10;
Elem.BendingAngle  	= A;
Elem.EntranceAngle 	= A1;
Elem.ExitAngle     	= A2;
Elem.ByError     	= 0;
Elem.K      			= K;

Elem.R1 = diag(ones(6,1));
Elem.R2 = diag(ones(6,1));
Elem.T1 = zeros(1,6);
Elem.T2 = zeros(1,6);

Elem.PolynomA			= [0 0 0 0];	 
Elem.PolynomB			= [0 K 0 0]; 

Elem.PassMethod 		= method;
