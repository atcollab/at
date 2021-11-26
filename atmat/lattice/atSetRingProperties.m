function ring = atSetRingProperties(ring,varargin)
%atSetRingProperties	Add or modify properties of the lattice
%
% newring=atSetRingProperties(ring [,key,value]...)
%   Add or modify the attributes of the RingParam element of the lattice,
%   Insert a new RingParam element if necessary
%
% Standard properties:
%   FamName         Name of the lattice
%   Energy          Ring energy in eV
%   Periodicity     Number of periods to make 2pi phase advance
%   Particle        particle (Particle object)
%   HarmNumber      Harmonic number
%
% Optional properties:
%   cavpts          Location of the main cavities (Used by atsetcavity)
%
% Additional custom fields may be added. They can de retrieved by
% atGetRingProperties and are saved in files.
%
% For fast access, the ring properties are stored in a RingParam element
% ideally located in the 1st position of the lattice. If there is no such
% element, atSetRingProperties will add it.
%
%See also atGetRingProperties

s=warning;                          % Save the warning state
warning('Off','AT:NoRingParam');    % Disable warning
[parmelem, idx] = atfindparam(ring, varargin{:});
warning(s);                         % Restore the warning state

if isempty(idx)
    % No RingParam element: insert a new one
    ring=[{parmelem};ring];
else
    % Update the existing RingParam element
    ring{idx}=parmelem;
end
end

