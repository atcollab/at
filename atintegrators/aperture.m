function z = aperture(fname,varargin)
%APERTURE Creates a aperture element in old a AT version (Obsolete)
%APERTURE('FAMILYNAME',limits, method)
%	creates a new family in the FAMLIST - a structure with fields
%		FamName			family name
%       Limits         (Xmin, Xmax, Ymin, Ymax) [m]
%		PassMethod		name of the function on disk to use for tracking
% returns assigned address in the FAMLIST that is uniquely identifies
% the family
%
%  NOTES
%  1. Obsolete: use ataperture instead
%  2. There is a better way to define aperture in each element
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector

ElemData=ataperture(fname,varargin{:});

global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;


