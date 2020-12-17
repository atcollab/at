function z=multipole(fname,L,PolynomA,PolynomB,method)
%MULTIPOLE Creates a thick multipole element in old AT versions (Obsolete)
% MULTIPOLE('FAMILYNAME',Length [m],PolynomA,PolynomB,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%	FamName			family name
%	Length			length[m]
%	ElemData.PolynomA= skew [dipole quad sext oct];	 
%	ElemData.PolynomB= normal [dipole quad sext oct]; 
%	PassMethod     name of the function to use for tracking
%
%   internally the additional structure fields are set:
%
%	NumIntSteps		Number of integration steps
%	MaxOrder
%	R1					6 x 6 rotation matrix at the entrance
%	R2        		6 x 6 rotation matrix at the entrance
%	T1					6 x 1 translation at entrance 
%	T2					6 x 1 translation at exit4
%
% returns assigned address in the FAMLIST that uniquely identifies
% the family
%
%  NOTES
%  1. Obsolete: use atmultipole or atthinmultipole instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

ElemData = atmultipole(fname,L,PolynomA,PolynomB,method);

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

