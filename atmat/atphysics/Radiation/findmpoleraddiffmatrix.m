function varargout=findmpoleraddiffmatrix(varargin) %#ok<STOUT>
%FINDMPOLERADDIFFMATRIX calculate radiation diffusion matrix of a multipole element
%DIFF=FINDMPOLERADDIFFMATRIX(ELEM,ORBIT_IN,ENERGY)
%
%ELEM:      AT multipole element
%ORBIT_IN:  input closed orbit
%ENERGY:    ring energy [eV]
%
%DIFF=FINDMPOLERADDIFFMATRIX(ELEM,ORBIT_IN)
%   takes energy from the 'Energy' field of the element
%
% for use in Ohmi's beam envelope formalism [1]
% See also OHMIENVELOPE
% [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
error('at:missingMex','missing MEX file.');
end
