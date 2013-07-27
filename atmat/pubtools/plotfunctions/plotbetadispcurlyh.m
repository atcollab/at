function plotdata=plotbetadispcurlyh(lindata,ring,dpp,varargin)

idx=cat(1,lindata.ElemIndex);
H=CurlyHlindata(lindata);

beta=cat(1,lindata.beta);                     % left axis
plotdata(1).values=beta;
plotdata(1).labels={'\beta_x','\beta_z'};
plotdata(1).axislabel='\beta [m]';
dispersion=cat(2,lindata.Dispersion)'; % right axis

plotdata(2).values=[dispersion(:,1)*100 H*10000];
plotdata(2).labels={'\eta_x','H*10^{2}'};
plotdata(2).axislabel='\eta_x, H [cm] ';
end
