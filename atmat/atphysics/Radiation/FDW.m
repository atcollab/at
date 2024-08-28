function varargout=FDW(varargin) %#ok<STOUT>
%FDW calculate radiation diffusion matrix of a wiggler element.
%DIFF=FDW(ELEM,ORBIT_IN,ENERGY)
%
%ELEM:      AT wiggler element
%ORBIT_IN:  input closed orbit
%ENERGY:    ring energy [GeV]
%
%DIFF=FDW(ELEM,ORBIT_IN)
%   takes energy from the 'Energy' field of the element
%
% for use in Ohmi's beam envelope formalism [1]
% See also OHMIENVELOPE
% [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
error('at:missingMex','missing MEX file.');
end