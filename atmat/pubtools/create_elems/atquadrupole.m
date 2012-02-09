function Elem=atquadrupole(fname,L,K,method)
%QUADRUPOLE('FAMILYNAME',Length [m],K,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields%		FamName    
%	FamName			family name
%	Length			length[m]
%	K				K-value of the quadrupole
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2        		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit
%	PassMethod     name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

Elem.FamName = fname;  % add check for existing identical family names
Elem.Length = L;
Elem.K         = K;
Elem.MaxOrder = 3;
Elem.NumIntSteps = 10;
Elem.PolynomA= [0 0 0 0];	 
Elem.PolynomB= [0 K 0 0];
Elem.R1 = diag(ones(6,1));
Elem.R2 = diag(ones(6,1));
Elem.T1 = zeros(1,6);
Elem.T2 = zeros(1,6);
Elem.PassMethod=method;