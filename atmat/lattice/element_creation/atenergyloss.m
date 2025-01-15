function elem = atenergyloss(fname, varargin)
%atenergyloss creates an energy loss element
%
%ELEM=ATENERGYLOSS(FAMNAME,ELOSS,PASSMETHOD)
%   FAMNAME:    family name
%   ELOSS:      Energy loss [eV]
%   PASSMETHOD: Tracking methods, defaults to 'IdentityPass'
%
%the "energy loss" element is taken into account in ATSUMMARY: it adds damping by
%contributing to the I2 integral, thus reducing the equilibrium emittance.
%But it does not generate any diffusion. This makes sense only if the losses
%summarised in the element occur in non-dispersive locations.

[rsrc,eloss,method]=decodeatargs({0,'IdentityPass'},varargin);
[eloss,rsrc]=getoption(rsrc,'EnergyLoss',eloss);
[method,rsrc]=getoption(rsrc,'PassMethod',method);
[cl,rsrc]=getoption(rsrc,'Class','EnergyLoss');
elem=atbaselem(fname,method,'Class',cl,'EnergyLoss',eloss,rsrc{:});
end
