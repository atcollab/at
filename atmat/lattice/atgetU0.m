function U0=atgetU0(ring)
%
% computes Energy loss per turn in eV .
%
%see also: ringpara atsetcavity

emass=510998.928/1e6; % electron mass in MeV

a = findcells(ring,'Energy');
if isempty(a);
   gamma = 3000/emass;
else
   gamma = ring{a(1)}.Energy/(emass*1e6); 
end

dpindex = findcells(ring,'BendingAngle');
lendp = getcellstruct(ring,'Length',dpindex); %bending magnet length
theta = getcellstruct(ring,'BendingAngle',dpindex); %bending angle
rho = lendp./theta;

I2 = sum(abs(theta./rho));

U0 = 14.085*(gamma*emass/1000)^4*I2*1000.; %eV

return
