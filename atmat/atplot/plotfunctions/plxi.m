function plotdata=plxi(lindata,ring,dpp) %#ok<INUSD>
%plotdata=plxi(lindata,ring,dpp)
%xi function in Touschek formula gives trans. velocity in beam frame
beta=cat(1,lindata.beta);
alpha=cat(1,lindata.alpha);
betax=beta(:,1);
betaz=beta(:,2);

delta=.03;

gamma=1.1820e+04;
epsx = 4e-9;

xi=delta^2/(gamma^2*epsx)*betax;

plotdata(1).values=xi;
plotdata(1).labels={'xi'};
plotdata(1).axislabel='xi';



plotdata(2).values=1./betax./betaz;
plotdata(2).labels={'1/(betax*betay)'};
plotdata(2).axislabel={'1/m^2'};

%plotdata(2).values=1./sigx./sigy;
%plotdata(2).labels={'1/(sigmax*sigmay)'};
%plotdata(2).axislabel={'1/m^2'};

end