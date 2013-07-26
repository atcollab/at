function plotdata=plBeamSize(lindata,r,dpp,varargin)
% plots H and V beam size
%
% created 21/09/2012

[r2,radindex,~]=atradon(r);
[envelope,~,~]=ohmienvelope(r2,radindex,1:length(r2)+1);
 
s=real(cat(2,envelope.Sigma))*1e3; % mum


%div=[[1;1], diff(s')'];

beta=cat(1,lindata.beta);                     % left axis
plotdata(1).values=beta;
plotdata(1).labels={'\beta_x','\beta_z'};
plotdata(1).axislabel='\beta [m]';

dispersion=cat(2,lindata.Dispersion)'; % right axis

plotdata(2).values=[dispersion(:,1) s'];% div'];
plotdata(2).labels={'\eta_x ','\sigma_x [m]','\sigma_y [m]'}; %,'\sigma_x'' [m]','\sigma_y'' [m]'};
plotdata(2).axislabel='dispersion [m] and \sigma [mm]';






return
