function [B66, M, rout] = findthinmpoleraddiffm(rin, PolynomA, PolynomB, L, irho, E0, max_order)
%FINDTHINMPOLERADDIFFM

% Physical constants used in calculations
persistent TWOPI CGAMMA M0C2 LAMBDABAR CER CU
if isempty(TWOPI) %Initialize constansts on the first call
    TWOPI       = 2*pi;
    CGAMMA      = 8.846056192e-05; 			% [m]/[GeV^3] Ref[1] (4.1)
    M0C2        = 5.10999060e5;              % Electron rest mass [eV]
    LAMBDABAR   = 3.86159323e-13;			% Compton wavelength/2pi [m]
    CER         = 2.81794092e-15;            % Classical electron radius [m]
    CU          = 1.323094366892892;         % 55/(24*sqrt(3))
end


% Calculate field from polynomial coefficients
P1 = i*PolynomA(1:max_order+1)+PolynomB(1:max_order+1);
Z1 = cumprod([1, (rin(1)+i*rin(3))*ones(1,max_order)]);
S1 = sum(P1.*Z1);
Bx = real(S1);
By = imag(S1);

B2P = B2perp([Bx; By+irho; 0], irho, rin);
B3P = B2P^(3/2);
p_norm = 1/(1+rin(5));
p_norm2 = p_norm^2;

CRAD = CGAMMA*E0^3/(TWOPI*1e27);
BB = CU * CER * LAMBDABAR *  (E0/M0C2)^5 * L * B3P * (1+rin(5))^4*...
				(1+rin(1)*irho + (rin(2)^2+rin(4)^2)*p_norm2/2);

% Propagate particle
rout = rin;

% Loss of energy (dp/p) due to radiation
rout(5) = rin(5) - CRAD*(1+rin(5))^2*B2P*...
    (1+rin(1)*irho + (rin(1)^2+rin(3)^2)*p_norm2/2)*L;

% Change in transverse momentum due to radiation
%   Angle does not change but dp/p changes due to radiation
%   and therefore transverse canonical momentum changes 
%   px = x'*(1+dp/p)
%   py = y'*(1+dp/p)
rout([2 4]) = rin([2 4])*(1+rout(5))/(1+rin(5));

% transverse kick due to magnetic field
rout(2) = rout(2) - L*(Bx-(rin(5)-rin(1)*irho)*irho);
rout(4) = rout(4) + L*By;

% pathlength
rout(6) = rout(6) + L*irho*rin(1); 


% Calculate transfer matrix at rin
P2 = i*PolynomA(2:max_order+1)+PolynomB(2:max_order+1);
Z2 = cumprod([1, (rin(1)+i*rin(3))*ones(1,max_order-1)]);
S2 = sum(P2.*(1:max_order).*Z2);

M = eye(6);
M(2,1)   = -L*real(S2);
M(2,3)   =  L*imag(S2);
M(4,1)   =  L*imag(S2);
M(4,3)   =  L*real(S2);
M(2,5)   =  L*irho;
M(2,1)   =  M(2,1) - L*irho*irho;
M(6,1)   =  L*irho;


%    Calculate Ohmi's diffusion matrix of a thin multipole  element 
%    For elements with straight coordinate system irho = 0
%    For curved elements the B polynomial (PolynomB in MATLAB) 
%    MUST NOT include the guide field  By0 = irho * E0 /(c*e)

B66 = zeros(6);
B66(2,2)    = BB*rin(2)^2*p_norm2;
B66(2,4)    = BB*rin(2)*rin(4)*p_norm2;
B66(4,2)    = B66(2,4);
B66(4,4)    = BB*rin(4)^2*p_norm2;
B66(5,2)    = BB*rin(2)*p_norm;
B66(2,5)    = B66(5,2);
B66(5,4)    = BB*rin(4)*p_norm;
B66(4,5)    = B66(5,4);
B66(5,5)    = BB;

function b2 = B2perp(B, irho, rin)
% Calculates sqr(|e x B|) , where 'e' is a unit vector in the direction of
% velocity. Components of the  velocity vector:
% ex = xpr; 
% ey = ypr; 
% ez = (1+x*irho);

E = [rin(2)/(1+rin(5));rin(4)/(1+rin(5));1+rin(1)*irho];
b2 = sum(cross(E/norm(E),B).^2);

