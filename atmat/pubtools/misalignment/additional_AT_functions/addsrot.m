function addsrot(elemindex, theta)
% ADDSROT increments the rotation obout the 's' axis (along the direction
% of particle motion) by an amount THETA in units of radians. ELEMINDEX is
% the index/indices of elements in THERING that are affected by this
% roataion. Clockwise rotation
%
% See also ADDXROT, ADDYROT

global THERING

numelems = length(elemindex);

if numelems ~= length(theta)
   error('ELEMINDEX and PSI must have the same number of elements');
end

for i = 1:numelems
   PolynomA = THERING{elemindex(i)}.PolynomA;
   PolynomB = THERING{elemindex(i)}.PolynomB;
   if isfield(THERING{elemindex(i)}, 'BendingAngle') && isfield(THERING{elemindex(i)}, 'EntranceAngle')
       PolynomB(0) = PolynomB(0) + THERING{elemindex(i)}.BendingAngle/THERING{elemindex(i)}.Length;
   end
   THERING{elemindex(i)}.PolynomA= imag(exp(complex(0,1)*theta(i))*complex(PolynomB,PolynomA));
   THERING{elemindex(i)}.PolynomB= real(exp(complex(0,1)*theta(i))*complex(PolynomB,PolynomA));
   if isfield(THERING{elemindex(i)}, 'BendingAngle') && isfield(THERING{elemindex(i)}, 'EntranceAngle')
       THERING{elemindex(i)}.PolynomB(0) = THERING{elemindex(i)}.PolynomB(0) - THERING{elemindex(i)}.BendingAngle/THERING{elemindex(i)}.Length;
   end
end