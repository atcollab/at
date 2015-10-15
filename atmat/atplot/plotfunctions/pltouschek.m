function plotdata=pltouschek(lindata,ring,dpp) %#ok<INUSD>
%plotdata=pltouschek(lindata,ring,dpp)
% plots curly H function and also 1/(sigx sigy) to understand local scattering rate
beta=cat(1,lindata.beta);
alpha=cat(1,lindata.alpha);
betax=beta(:,1);
betaz=beta(:,2);
alphax=alpha(:,1);
gammax=(1+alphax.*alphax)./betax;
eta=cat(2,lindata.Dispersion)';
etax=eta(:,1);
etaxp=eta(:,2);
Hx=gammax.*etax.*etax+2*alphax.*etax.*etaxp+betax.*etaxp.*etaxp;

%now compute beamenvelopes
%atxdata=atx(ring,dpp,1:length(ring)+1);
%sig=cat(3,atxdata.beam66);
%sigx=sig(1,1,:);
%sigx=squeeze(sigx);
%sigy=sig(3,3,:);
%sigy=squeeze(sigy);

%momap_h = varargin(1);

plotdata(1).values=Hx;
plotdata(1).labels={'H_x'};
plotdata(1).axislabel='H [m]';

plotdata(2).values=1./betax./betaz;
plotdata(2).labels={'1/(betax*betay)'};
plotdata(2).axislabel={'1/m^2'};

%plotdata(2).values=1./sigx./sigy;
%plotdata(2).labels={'1/(sigmax*sigmay)'};
%plotdata(2).axislabel={'1/m^2'};

end