function plotdata=plotsqrtbetadispcurlyh(lindata,ring,dpp,varargin)

H=CurlyHlindata(lindata);

beta=cat(1,lindata.beta);                     % left axis
plotdata(1).values=sqrt(beta);
plotdata(1).labels={'\beta_x^{0.5}','\beta_z^{0.5}'};
plotdata(1).axislabel='\beta^{0.5} [m^{0.5}]';
dispersion=cat(2,lindata.Dispersion)'; % right axis

plotdata(2).values=[dispersion(:,1)*100 H*10000];
plotdata(2).labels={'\eta_x [cm]','H [10-4]'};
plotdata(2).axislabel='dispersion [cm]';
end
