function [I1,I2,I3,I4,I5] = WigglerRadiation(ring,lindata)
%WIGGLERRADIATION	Compute the radiation integrals in wigglers
%
%[I1,I2,I3,I4,I5] = WIGGLERRADIATION(RING,LINDATA)
%
%RING       Lattice structure
%LINDATA    Output of atlinopt for all lattice elements
%
%WigglerRadiation computes the radiation integrals for all wigglers with
%the following approximations:
%
%- The self-induced dispersion is neglected in I4 and I5, but is is used as
%  a lower limit for the I5 contribution
%
% I1, I2 are integrated analytically
% I3 is integrated analytically for a single harmonic, numerically otherwise

nstep=60;
e_mass=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;	% eV
cspeed = PhysConstant.speed_of_light_in_vacuum.value;                   % m/s

iswiggler=atgetcells(ring,'Class','Wiggler');
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
    [I1,I2,I3,I4,I5]=deal(0);
end

    function [di1,di2,di3,di4,di5]=wigrad(elem,dini)
        le=elem.Length;
        alphax0=dini.alpha(1);
        betax0=dini.beta(1);
        gammax0=(alphax0.*alphax0+1)./betax0;
        eta0 = dini.Dispersion(1);
        etap0 = dini.Dispersion(2);
        H0=gammax0*eta0*eta0 + 2*alphax0*eta0*etap0 + betax0*etap0*etap0;
        avebetax=betax0+alphax0*le+gammax0*le*le/3;

        kw=2*pi/elem.Lw;
        rhoinv=elem.Bmax/Brho;
        coefh=elem.By(2,:);
        coefv=elem.Bx(2,:);
        coef2=[coefh coefv];
        if length(coef2)==1     % Analytical I3
            di3=(rhoinv*coef2)^3*4*le/3/pi;
        else                    % Numerical I3
            [bx,bz]=Baxis(elem,linspace(0,elem.Lw,nstep+1));
            B2=bx.*bx+bz.*bz;
            rinv=sqrt(B2)/Brho;
            di3=trapz(rinv.^3)*le/nstep;
        end
        di2=elem.Length*(coefh*coefh'+coefv*coefv')*rhoinv*rhoinv/2;
        di1=-di2/kw/kw;
        di4=0;
        d5lim=4*avebetax*le*rhoinv^5/15/pi/kw/kw;
        di5=max(H0*di3,d5lim);
%         fprintf('%s\t%e\t%e\t%e\t%e\t%e\t(%e,%e,%e)\n',elem.FamName,di1,di2,di3,di4,di5,d5lim,H0*di3,dini.Dispersion(1));
    end

    function [bx,bz] = Baxis(wig,s)
        %Bwig   Compute the field on the axis of a generic wiggler
        %
        %The field of an horizontal wiggler is represented by a sum of harmonics,
        %each being described by a 6x1 column in a matrix By such as:
        %
        %   Bz/Bmax = -By2 * cos(By5*kw*s + By6)
        %
        %The field of an vertical wiggler is represented by a sum of harmonics,
        %each being described by a 6x1 column in a matrix Bx such as:
        %
        %   Bx/Bmax =  Bx2 * cos(Bx5*kw*s + Bx6)
        
        kw=2*pi/wig.Lw;
        Bmax=wig.Bmax;
        kws=kw*s;
        [bxh,bzh]=cellfun(@harmh,num2cell(wig.By,1),'UniformOutput',false);
        [bxv,bzv]=cellfun(@harmv,num2cell(wig.Bx,1),'UniformOutput',false);
        bx=sum(cat(3,bxh{:},bxv{:}),3);
        bz=sum(cat(3,bzh{:},bzv{:}),3);
        
        function [bx,bz]=harm(pb)
            bz=-Bmax*pb(2)*cos(pb(5)*kws+pb(6));
            bx=zeros(size(kws));
        end
        function [bx,bz]=harmh(pb)
            [bx,bz]=harm(pb);
        end
        function [bx,bz]=harmv(pb)
            [bz,fz]=harm(pb);
            bx=-fz;
        end
    end

end

