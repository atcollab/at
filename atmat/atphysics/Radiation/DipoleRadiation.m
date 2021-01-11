function [I1,I2,I3,I4,I5,I6,Iv] = DipoleRadiation(ring,lindata)
%DIPOLERADIATION	Compute the radiation integrals in dipoles

% Analytical integration from:
%
% EVALUATION OF SYNCHROTRON RADIATION INTEGRALS
% R.H. Helm, M.J. Lee, P.L. Morton and M. Sands
% SLAC-PUB-1193, March 1973


angle=atgetfieldvalues(ring,'BendingAngle');
isdipole=isfinite(angle) & (angle~=0);
vini=lindata([isdipole;false])';
vend=lindata([false;isdipole])';
[di1,di2,di3,di4,di5,di6,div]=cellfun(@diprad,ring(isdipole),num2cell(vini),num2cell(vend));
I1=sum(di1);
I2=sum(di2);
I3=sum(di3);
I4=sum(di4);
I5=sum(di5);
I6=sum(di6);
Iv=sum(div);

    function [di1,di2,di3,di4,di5,di6,div]=diprad(elem,dini,dend)
            betax0 = dini.beta(1);
            betaz0 = dini.beta(2);
            alphax0 = dini.alpha(1);
            eta0 = dini.Dispersion(1);
            etap0 = dini.Dispersion(2);

            ll = elem.Length;
            theta = elem.BendingAngle;
            rho = ll / theta;
            rho2 = rho * rho;
            K = getfoc(elem);
            kx2 = K + 1.0/rho2;
            kz2 = -K;
            eps1 = tan(elem.EntranceAngle) / rho;
            eps2 = tan(elem.ExitAngle) / rho;

            eta3 = dend.Dispersion(1);
            alphax1 = alphax0 - betax0*eps1;
            gammax1 = (1.0 + alphax1 * alphax1) / betax0;
            alphaz1 = dini.alpha(2) + betaz0*eps1;
            alphaz2 = dend.alpha(2) - dend.beta(2)*eps2;
            gammaz1 = (1.0 + alphaz1*alphaz1) / betaz0;
            etap1 = etap0 + eta0*eps1;
            etap2 = dend.Dispersion(2) - eta3*eps2;

            h0 = gammax1*eta0*eta0 + 2.0*alphax1*eta0*etap1 + betax0*etap1*etap1;

            if kx2 ~= 0.0
                if kx2 > 0.0         % Focusing
                    kl = ll * sqrt(kx2);
                    ss = sin(kl) / kl;
                    cc = cos(kl);
                else                % Defocusing
                    kl = ll * sqrt(-kx2);
                    ss = sinh(kl) / kl;
                    cc = cosh(kl);
                end
                eta_ave = (theta - (etap2 - etap1)) / kx2 / ll;
                bb = 2.0 * (alphax1*eta0 + betax0*etap1) * rho;
                aa = -2.0 * (alphax1*etap1 + gammax1*eta0) * rho;
                h_ave = h0 + (aa * (1.0-ss) + bb * (1.0-cc) / ll ...
                              + gammax1 * (3.0-4.0*ss+ss*cc) / 2.0 / kx2 ...
                              - alphax1 * (1.0-cc)^2 / kx2 / ll ...
                              + betax0 * (1.0-ss*cc) / 2.0 ...
                              ) / kx2 / rho2;
            else
                eta_ave = 0.5 * (eta0 + eta3) - ll*ll / 12.0 / rho;
                hp0 = 2.0 * (alphax1 * eta0 + betax0 * etap1) / rho;
                h2p0 = 2.0 * (-alphax1*etap1 + betax0/rho - gammax1*eta0) / rho;
                h_ave = h0 + hp0*ll/2.0 + h2p0*ll*ll/6.0 ...
                    - alphax1*ll^3/4.0/rho2 ...
                    + gammax1*ll^4/20.0/rho2;
            end
            if kz2 ~= 0.0
                bz_ave=(gammaz1 + kz2*betaz0 + (alphaz2-alphaz1)/ll)/2/kz2;
            else
                bz_ave=betaz0-alphaz1*ll + gammaz1*ll*ll/3;
            end

            di1 = eta_ave * ll / rho;
            di2 = ll / rho2;
            di3 = ll / abs(rho) / rho2;
            di4 = eta_ave * ll * (2.0*K+1.0/rho2) / rho ...
                - (eta0*eps1 + eta3*eps2)/rho;
            di5 = h_ave * ll / abs(rho) / rho2;
            di6 = kz2^2*eta_ave^2*ll;
            div = abs(theta)/rho2*bz_ave;
    end

    function K=getfoc(elem)
        if isfield(elem,'PolynomB') && length(elem.PolynomB) >= 2
            K = elem.PolynomB(2);
            if isfield(elem,'K') && elem.K ~= K
                warning('AT:InconsistentK',...
                'Values in K and PolynomB(2) are different. Using PolynomB(2)');
            end
        elseif isfield(elem,'K')
            K = elem.K;
        else
            K = 0;
        end
    end
end

