function U0 = atgetU0(ring)
%ATGETU0 - computes Energy loss per turn in eV .
%
%  INPUTS 
%    1. ring ring structure
%
%  OUTPUTS
%    1. U0 energy loss per turn in eV
%
%
% See also  ringpara atsetcavity atenergy

%
%% S.Liuzzo modified to use atenergy.m

[~,~,~,~,U0]=atenergy(ring);

return

% 
% emass=510998.928/1e6; % electron mass in MeV
% 
% a = findcells(ring,'Energy');
% if isempty(a);
%    gamma = 3000/emass;
% else
%    gamma = ring{a(1)}.Energy/(emass*1e6); 
% end
% 
% dipindex = find(atgetcells(ring,'BendingAngle'));
% lendp = atgetfieldvalues(ring,dipindex,'Length'); %bending magnet length
% theta = atgetfieldvalues(ring,dipindex,'BendingAngle'); %bending angle
% rho = lendp./theta;
% 
% I2 = sum(abs(theta./rho));
% 
% U0 = 14.085*(gamma*emass/1000)^4*I2*1000.; %eV
% 
% 
% 
% 
% return
