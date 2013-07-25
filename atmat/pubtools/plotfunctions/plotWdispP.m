function plotdata=plotWdispP(lindata,ring,dpp,varargin)

idx=cat(1,lindata.ElemIndex);
DE=0.001;

[lz,~,~]=atlinopt(ring,0,idx);
[lpd,~,~]=atlinopt(ring,DE,idx);
[lmd,~,~]=atlinopt(ring,-DE,idx);
bz=cat(1,lz.beta);
bp=cat(1,lpd.beta);
bm=cat(1,lmd.beta);   % left axis
az=cat(1,lz.alpha);
ap=cat(1,lpd.alpha);
am=cat(1,lmd.alpha);   % left axis

aa=(bp-bm)./2./DE./bz;
bb=(ap-am)./2./DE-az./bz.*(bp-bm)./2./DE;
size(aa)
size(bb)
plotdata(1).values=sqrt(aa.^2+bb.^2);

plotdata(1).labels={'\beta_x/D\delta','\beta_z/D\delta'};
plotdata(1).axislabel='Wx Wy [m]';

dp=cat(2,lpd.Dispersion)'; % right axis
dm=cat(2,lmd.Dispersion)';
plotdata(2).values=[(dp-dm)./2/DE*100];
plotdata(2).labels={'D\eta_{x}/D\delta [cm]'};
plotdata(2).axislabel='D'' [cm]';
end
