function [props, idx] = atGetRingProperties(ring)
%ATGETRINGPROPERTIES Get the ring properties
%
% properties=ATGETRINGPROPERTIES(ring)
%	Extract data from the RingParam element of the lattice, if present,
%	or from the lattice elements
%
% ring:             Ring structure
%
% properties: Structure with fields:
%   FamName         Name of the lattice
%   Energy          Ring energy in eV
%   Periodicity     Number of periods to make 2pi phase advance
%   Particle        particle (Particle object)
%   HarmNumber      Harmonic number
%
% Optional fields:
%   cavpts          Location of the main cavities (Used by atsetcavity)
%
% For fast access, the ring properties are stored in a RingParam element
% ideally located in the 1st position of the lattice. Without such element,
% the properties are deduced from the lattice contents. This is much slower
% and ATGETRINGPROPERTIES displays a warning indicating how to add the
% RingParam element:
%
%>>ring=atSetRingProperties(ring)
%
%  See also ATSETRINGPROPERTIES

[parmelem, idx] = atfindparam(ring);
particle = atparticle.loadobj(parmelem.Particle);
props=rmfield(parmelem,{'Length','Class','PassMethod','Particle'});
props.Particle = particle;
end
