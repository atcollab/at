function Elem=atsextupole(fname,L,S,method)

%SEXTUPOLE('FAMILYNAME',Length [m],S,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields%		FamName    
%	FamName			family name
%	Length			length[m]
%	S					S-strngth of the sextupole
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2        		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit4
%	ElemData.PolynomA= [0 0 0 0];	 
%	ElemData.PolynomB= [0 0 S 0]; 
%	PassMethod     name of the function to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family


Elem.FamName = fname;  % add check for identical family names
Elem.Length = L;
Elem.MaxOrder = 3;
Elem.NumIntSteps = 10;
Elem.R1 = diag(ones(6,1));
Elem.R2 = diag(ones(6,1));
Elem.T1 = zeros(1,6);
Elem.T2 = zeros(1,6);
Elem.PolynomA= [0 0 0 0];	 
Elem.PolynomB= [0 0 S 0]; 
Elem.PassMethod=method;
