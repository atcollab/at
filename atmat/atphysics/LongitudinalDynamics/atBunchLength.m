function BL = atBunchLength (ring,Ib,Zn)
% bunch length due to the potential well effect
% the output is the zerocurrent bunch length x bunch lengthening
%
%   BL = atBunchLength(ring,Ib,Zn)
%
% Ib is the bunch current [A] (it may be a vector for multiple values)
% Zn is the longitudinal broadband impedance [Ohms]
% ring is the at ring without radiation
% BL is the bunch length in metres 
%
%   see also: BunchLength

rp=ringpara(ring);
[E0,~,Vrf,h,U0]=atenergy(ring);
alpha=rp.alphac;
sigdelta=rp.sigma_E;
circ=findspos(ring,length(ring)+1);

BL = BunchLength(Ib,Zn,Vrf,U0,E0,h,alpha,sigdelta,circ);

end
