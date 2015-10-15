function plotdata=plEmitContrib(lindata,ring,dpp,varargin)
% plots H/rho³ at every dipole

idx=cat(1,lindata.ElemIndex);
H=CurlyH(ring,dpp,idx(:)');

r=zeros(size(H));
bend=findcells(ring,'BendingAngle');
r(bend)=getcellstruct(ring,'Length',bend)./getcellstruct(ring,'BendingAngle',bend);
emitcontr=H./(r.^3)*1e9;
emitcontr(isinf(emitcontr))=0;

beta=cat(1,lindata.beta);                     % left axis
plotdata(1).values=beta;
plotdata(1).labels={'\beta_x','\beta_z'};
plotdata(1).axislabel='\beta [m]';

dispersion=cat(2,lindata.Dispersion)'; % right axis
plotdata(2).values=[dispersion(:,1)*100 emitcontr H*10000];
plotdata(2).labels={'\eta_x cm','H/r^{3}*1e9','H 10-4'};
plotdata(2).axislabel='dispersion [cm]';

end
