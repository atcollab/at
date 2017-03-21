function plotdata=plotERAperture(lindata,ring,dpp,varargin) %#ok<INUSD>
%function plotdata=plotERAperture(lindata,ring,dpp,varargin) %#ok<INUSD>
%
% Plots RApertures EApertures aperture.
% 
% use with atplot(atlattice,@plotAperture);
%

rapind=findcells(ring,'RApertures');
eapind=findcells(ring,'EApertures');

xm=getcellstruct(ring,'RApertures',rapind,1);
ym=getcellstruct(ring,'RApertures',rapind,3);
xp=getcellstruct(ring,'RApertures',rapind,2);
yp=getcellstruct(ring,'RApertures',rapind,4);
eh=getcellstruct(ring,'EApertures',eapind,1);
ev=getcellstruct(ring,'EApertures',eapind,2);

Xp=[nan(size(ring)); nan];
Xm=Xp;
Yp=Xp;
Ym=Xp;
Eh=Xp;
Ev=Xp;

Xp(rapind)=xp;
Xm(rapind)=xm;
Yp(rapind)=yp;
Ym(rapind)=ym;
Eh(eapind)=eh;
Ev(eapind)=ev;

Xp=fixgaps(Xp);
Xm=fixgaps(Xm);
Yp=fixgaps(Yp);
Ym=fixgaps(Ym);
Eh=fixgaps(Eh);
Ev=fixgaps(Ev);

plotdata(1).values=[Xp Xm Yp Ym]*1e2;%
plotdata(1).labels={'x ','-x','y','-y'};
plotdata(1).axislabel='rectangular aperture [cm]';
% 
plotdata(2).values=[Eh Ev]*1e2;%
plotdata(2).labels={'hor.','ver.'};
plotdata(2).axislabel='elliptic aperture [cm]';
% 

end 

function y=fixgaps(x)
% FIXGAPS Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN
% in the input time series (may be complex), but ignores
% trailing and leading NaN.
%

% R. Pawlowicz 6/Nov/99

y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0;


y(bd)=interp1(gd,x(gd),find(bd)); 

end