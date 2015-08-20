function x = atbeam(np,sigma,orbit)
%PARTICLES=ATBEAM(NP,SIGMA)  Generate a particle distribution according to a sigma matrix
%
%NP:    number of particles
%SIGMA: beam matrix (2x2, 4x4, 6x6)
%
%PARTICLES=ATBEAM(NP,SIGMA,ORBIT) adds a center of mass to the distribution
%
%See also: ATPLOTBEAM

% ampl=sqrt(-2*log(rand(3,np)));
% phase=2*pi()*rand(3,np);
% v=[ampl.*cos(phase);ampl.*sin(phase)];
v=randn(size(sigma,1),np);
try
    l=chol(sigma);
catch
    a=[chol(sigma([1 2 5 6],[1 2 5 6])) zeros(4,2);zeros(2,6)];
    l=a([1 2 5 6 3 4],[1 2 5 6 3 4]);
end
if nargin < 3
    x=l'*v;
else
    x=orbit(:,ones(1,np)) + l'*v;
end
end
