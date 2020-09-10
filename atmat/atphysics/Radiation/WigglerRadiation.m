function [I1,I2,I3,I4,I5] = WigglerRadiation(ring,lindata)
%DIPOLERADIATION	Compute the radiation integrals in dipoles

e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s

iswiggler=atgetfieldvalues(ring,'Class','Wiggler');
if any(iswiggler)
    energy=unique(atgetfieldvalues(ring(iswiggler),'Energy'));
    if length(energy) > 1
        error('AT:NoEnergy','Energy field not equal for all elements');
    end
    gm=energy/e_mass;
    bt=sqrt(1.0-1.0/gm/gm);
    Brho=bt*energy/cspeed;
    vini=lindata([iswiggler;false])';
    [di1,di2,di3,di4,di5]=cellfun(@wigrad,ring(iswiggler),num2cell(vini));
    I1=sum(di1);
    I2=sum(di2);
    I3=sum(di3);
    I4=sum(di4);
    I5=sum(di5);
else
    I1=0;
    I2=0;
    I3=0;
    I4=0;
    I5=0;
end

    function [di1,di2,di3,di4,di5]=wigrad(elem,dini)
        rhoinv=elem.Bmax/Brho;
        coefh=elem.By(2,:);
        coefv=elem.Bx(2,:);
        [bx,bz,~]=Bwig(elem,0,0,linspace(0,elem.Lw,61));
        B2=bx.*bx+bz.*bz;
        rinv=sqrt(B2)/Brho;
        di2=trapz(B2)*elem.Length/Brho;
        di3=trapz(rinv.^3)*elem.Length;
        di2=elem.Length*(coefh*coefh'+coefv*coefv')*rhoinv*rhoinv/2;
        betax0 = dini.beta(1);
        alphax0 = dini.alpha(1);
        gammax0 = (1+alphax0*alphax0)/betax0;
        eta0 = dini.Dispersion(1);
        etap0 = dini.Dispersion(2);
        H0=gammax0*eta0*eta0 + 2*alphax0*eta0*etap0 + betax0*etap0*etap0;
        di1=0;
        di4=0;
        di5=H0*di2;
    end
end

