function plotdata=plotB0curlyh(lindata,ring,dpp,varargin)

spl=299792458;
E0 = atenergy(ring);
Brho=E0/spl;

idx=cat(1,lindata.ElemIndex);
H=CurlyHlindata(lindata);
B0=zeros(size(H));
dips=findcells(ring,'BendingAngle');
B0(dips)=Brho*getcellstruct(ring,'BendingAngle',dips)./...
    getcellstruct(ring,'Length',dips);


beta=cat(1,lindata.beta);                     % left axis
plotdata(1).values=B0;
plotdata(1).labels={['B [T] ' num2str(E0*1e-9) 'GeV']};
plotdata(1).axislabel='B [T]';
dispersion=cat(2,lindata.Dispersion)'; % right axis

plotdata(2).values=[dispersion(:,1)*100 H*10000];
plotdata(2).labels={'\eta_x [cm]','H [10-2]'};
plotdata(2).axislabel='dispersion , H [cm]';
end
