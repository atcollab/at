function plotdata=plPolynomBSxtOct(lindata,ring,~,varargin)
%PLPOLYNOMBSXTOCT Plots Bn for sextupole and octupole magnets
%DEFAULTPLOT    Default plotting function for ATPLOT
%
%Plots polynomB for ring and ring1

CoD=cat(2,lindata.ClosedOrbit);
% PolyBVal=zeros(size(CoD(1,:)));
% PolyBVal1=zeros(size(CoD(1,:)));

% PolynomBVal1=zeros(size(CoD(1,:)));
% PolynomBVal2=zeros(size(CoD(1,:)));
PolynomBVal3=zeros(size(CoD(1,:)));
PolynomBVal4=zeros(size(CoD(1,:)));
ind=findcells(ring,'PolynomB');

% PolynomBVal1(ind)=getcellstruct(ring,'PolynomB',ind,1,1);
% PolynomBVal2(ind)=getcellstruct(ring,'PolynomB',ind,1,2);
PolynomBVal3(ind)=getcellstruct(ring,'PolynomB',ind,1,3);
PolynomBVal4(ind)=getcellstruct(ring,'PolynomB',ind,1,4);


plotdata(1).values=  PolynomBVal3';
plotdata(2).values=  PolynomBVal4';

plotdata(1).labels={'b_{3}'};
plotdata(1).axislabel='b_{3} [1/m^{2}]';
plotdata(2).labels={'b_{4}'};
plotdata(2).axislabel='b_{4} [1/m^{3}]';
% dispersion=cat(2,lindata.Dispersion)';
% plotdata(2).values=dispersion(:,3);
% plotdata(2).labels={'\eta_y'};
% plotdata(2).axislabel='vertical dispersion [m]';

end
