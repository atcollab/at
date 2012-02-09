function Elem=solenoid(fname,L,KS,method)

% z=solenoid('FAMILYNAME',Length [m],KS,'METHOD')
%	creates a new family in the FAMLIST - a structure with field
%	FamName			family name
%	Length			length[m]
%	KS              solenoid strength KS [rad/m]
%	PassMethod     name of the function to use for tracking
%
%   function returns assigned address in the FAMLIST that uniquely identifies
%   the family
%
%   Additional structures being set up (initialized to default values within this routine):   
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2           		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit

Elem.FamName = fname;  % add check for existing identical family names
Elem.Length = L;
Elem.K         = KS;
Elem.MaxOrder = 3;
Elem.NumIntSteps = 10;
Elem.R1 = diag(ones(6,1));
Elem.R2 = diag(ones(6,1));
Elem.T1 = zeros(1,6);
Elem.T2 = zeros(1,6);
Elem.PassMethod=method;