function z = drift(fname,L,method)
%DRIFT Creates a drift element in old a AT version (Obsolete)
%DRIFT('FAMILYNAME',Length [m],'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length			length[m]
%		PassMethod		name of the function on disk to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
%
%  NOTES
%  1. Obsolete: use atdrift instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

ElemData = atdrift(fname,L,method);

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;


