function [Bcum, Mcum, r] = findthickmpoleraddifm(rin, PolynomA, PolynomB,L, irho, E0, max_order,num_steps)
%FINDTHICKMPOLERADDIFFM

% Fourth order-symplectic integrator constants
persistent DRIFT1 DRIFT2 KICK1 KICK2
if isempty(DRIFT1)
    DRIFT1   = 0.6756035959798286638;
    DRIFT2   = -0.1756035959798286639;
    KICK1    =  1.351207191959657328;
    KICK2    = -1.702414383919314656;
end



SL = L/num_steps;
L1 = SL*DRIFT1;
L2 = SL*DRIFT2;
K1 = SL*KICK1;
K2 = SL*KICK2;

Mcum = eye(6);
Bcum = zeros(6);
r = rin;

for m=1:num_steps % Loop over slices
    
    [M, r] = driftm66(L1,r);
    Bcum = M*Bcum*M';
    Mcum = M*Mcum;
    
    [B, M, r] = findthinmpoleraddiffm(r, PolynomA, PolynomB, K1, irho, E0, max_order);
    Bcum = M*Bcum*M' + B;
    Mcum = M*Mcum;
    
    [M, r] = driftm66(L2,r);
    Bcum = M*Bcum*M';
    Mcum = M*Mcum;
    
    [B, M, r] = findthinmpoleraddiffm(r, PolynomA, PolynomB, K2, irho, E0, max_order);
    Bcum = M*Bcum*M' + B;
    Mcum = M*Mcum;
				
    [M, r] = driftm66(L2,r);
    Bcum = M*Bcum*M';
    Mcum = M*Mcum;
	
    [B, M, r] = findthinmpoleraddiffm(r, PolynomA, PolynomB, K1, irho, E0, max_order);
    Bcum = M*Bcum*M' + B;
    Mcum = M*Mcum;
    
    [M, r] = driftm66(L1,r);
    Bcum = M*Bcum*M';
    Mcum = M*Mcum;

end


function [M, rout] = driftm66(L,r);
% transfer matrix of a drift - map linearized at r

Pnorm = 1/(1+r(5));
NormL = L*Pnorm;
M = eye(6);
M([7 21]) = NormL;


M([1 3],5) = -NormL*r([2,4])*Pnorm;
M(6,[2 4]) = -M([1 3],5)';
M(6,5) =  -NormL*Pnorm*sum(r([2,4]).^2);

rout = r;
rout([1 3]) = r([1 3]) + r([2 4])*NormL;
rout(6) = r(6) + NormL*Pnorm*sum(r([2,4]).^2)/2;

