function setyrot(elemindex, theta)
% ADDSROT increments the rotation obout the 's' axis (along the direction
% of particle motion) by an amount THETA in units of radians. ELEMINDEX is
% the index/indices of elements in THERING that are affected by this
% roataion. Clockwise rotation
%
% See also ADDSROT, ADDXROT, SETSROT, SETXROT, SETSHIFT, ADDSHIFT, ADDYROT

global THERING

numelems = length(elemindex);

if numelems ~= length(theta)
    error('ELEMINDEX and THETA must have the same number of elements');
end

T = tan(theta);

for i = 1:numelems
    elemlen = THERING{elemindex(i)}.Length;
    % Add x translation and kick to particle coordiante to simulate
    % rotation of element about y-axis.
    THERING{elemindex(i)}.T1(1) = elemlen*0.5*T(i);
    THERING{elemindex(i)}.T1(2) = -theta(i);
    THERING{elemindex(i)}.T2(1) = -elemlen*0.5*T(i);
    THERING{elemindex(i)}.T2(2) = theta(i);
end