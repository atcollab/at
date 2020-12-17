function z = corrector(fname,L,kickangle, method)
%CORRECTOR Creates a corrector element in old a AT version (Obsolete)
%CORRECTOR('FAMILYNAME',LENGTH,ANGLE,'METHOD')
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%		Length 			is set to 0 for  marker type 
%		KickAngle       [kickx, kicky] in radians (small) - unis of d(x,y)/ds
%       PassMethod		name of the function on disk to use for tracking
%
% returns assigned index in the FAMLIST that is uniquely identifies
% the family
%
%  NOTES
%  1. Obsolete: use atcorrector instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

ElemData = atcorrector(fname,L,kickangle,method);

global FAMLIST
z = length(FAMLIST)+1; % number of declared families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;

