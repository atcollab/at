function [I1,I2,I3,I4,I5] = ElossRadiation(ring,~)
%ELOSSRADIATION Compute the radiation integrals in EnergyLoss elements
%
%[I1,I2,I3,I4,I5] = ELOSSRADIATION(RING,LINDATA)
%
%RING       Lattice structure
%LINDATA    Output of atlinopt for all lattice elements (not used)

% Losses = Cgamma/2/pi*EGeV^4*I2
energy=atGetRingProperties(ring,'Energy');
e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
e_radius=PhysConstant.classical_electron_radius.value;                  % m
Cgamma=4.0*pi*e_radius/3/e_mass^3;                                      % [m/eV^3]
coef=Cgamma/2/pi*energy^4;

[I1,I3,I4,I5]=deal(0);
eloss=atgetcells(ring,'Class','EnergyLoss');
di2=cellfun(@eloss_i2,ring(eloss));
I2=sum(di2);
    
    function i2=eloss_i2(elem)
        i2=elem.EnergyLoss/coef;
    end

end