function plotdata=plPolynomBComp(lindata,ring,dpp,varargin) %#ok<INUSD>
% PLPOLYNOMBCOMP default plotting function for ATPLOT
%
% Plots polynomB for ring and ring1

CoD=cat(2,lindata.ClosedOrbit);
PolyBVal=zeros(size(CoD(1,:)));
PolyBVal1=zeros(size(CoD(1,:)));

PolynomBVal1=zeros(size(CoD(1,:)));
PolynomBVal2=zeros(size(CoD(1,:)));
PolynomBVal3=zeros(size(CoD(1,:)));
PolynomBVal4=zeros(size(CoD(1,:)));
ind=findcells(ring,'PolynomB');

PolynomBVal1(ind)=getcellstruct(ring,'PolynomB',ind,1,1);
PolynomBVal2(ind)=getcellstruct(ring,'PolynomB',ind,1,2);
PolynomBVal3(ind)=getcellstruct(ring,'PolynomB',ind,1,3);
PolynomBVal4(ind)=getcellstruct(ring,'PolynomB',ind,1,4);


plotdata(1).values=[PolynomBVal1' PolynomBVal2'/10 ...];%...
                     PolynomBVal3'/100 PolynomBVal4'/5000 ];
plotdata(1).labels={'dipole','quadrupole/10',...
                    'sextuople/100','octupole/5000'};
plotdata(1).axislabel='PolynomB';
% dispersion=cat(2,lindata.Dispersion)';
% plotdata(2).values=dispersion(:,3);
% plotdata(2).labels={'\eta_y'};
% plotdata(2).axislabel='vertical dispersion [m]';

end
