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
%   see also: BunchLength, blgrowth

is6d=check_6d(ring);
circ=findspos(ring,length(ring)+1);

if ~is6d
    [E0,~,Vrf,h,U0]=atenergy(ring);
    rp=ringpara(ring);
    alpha=rp.alphac;
    sigdelta=rp.sigma_E;
    BL = BunchLength(Ib,Zn,Vrf,U0,E0,h,alpha,sigdelta,circ);
else
    [~,ringdata]=atx(ring,1);
    [~,~,Vrf,h]=atenergy(ring);
    alpha = ringdata.alpha;
    E0 = ringdata.energy;
    U0 = ringdata.eloss;
    sigdelta = ringdata.espread;
    bl0 = ringdata.blength;
    BL = bl0 * blgrowth(Ib,Zn,Vrf,U0,E0,h,alpha,sigdelta);
end
end
