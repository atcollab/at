function [z] = wiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method,varargin)
% wiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method,energy)
%
% FamName	family name
% Ltot		total length of the wiggle
% Lw		total length of the wiggle
% Bmax	 	peak wiggler field [Tesla]
% Nstep		num of integration steps per period
% Nmeth		symplectic integration method, 2nd or 4th order: 2 or 4
% By		wiggler harmonics for horizontal wigglers
% Bx		wiggler harmonics for vertical wigglers
% method        name of the function to use for tracking
% energy    Beam Energy: Default is GLOBVAL.E0 in case it exists 0 otherwise
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

%  NOTES
%  1. Obsolete: use atwiggler instead
%
%  See also atdrift, atquadrupole, atsextupole, atsbend, atskewquad,
%          atmultipole, atthinmultipole, atmarker, atcorrector
global FAMLIST GLOBVAL

energy=0;
if numel(varargin)==0
    if isfield(GLOBVAL,'E0')
        energy=GLOBVAL.E0;
    end
elseif numel(varargin)==1
    energy=varargin{1};
end

ElemData = atwiggler(fname, Ltot, Lw, Bmax, energy, method,'Nstep',Nstep, 'Nmeth',Nmeth, 'By',By, 'Bx', Bx);


z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;
