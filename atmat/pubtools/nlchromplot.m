function nlchromplot(ring,min_delta,max_delta,steps,nper,fig,dash)
%example:   nlchromplot(esrf,-.04,.04,30,16,1)
delta_dif=(max_delta-min_delta)/steps;
delta_vals=linspace(min_delta,max_delta,steps);
nu_delta=zeros(2,steps);
for j=1:steps
[lindata,nu_delta(:,j)]=atlinopt(ring,min_delta+(j-1)*delta_dif);
%[lindata,nu_delta(:,j)]=linopt(ring,min_delta+(j-1)*delta_dif);
end
nux=nper*nu_delta(1,:);
nux=nux-floor(nux);
nuz=nper*nu_delta(2,:);
%nuz=1-(nuz-floor(nuz));
nuz=nuz-floor(nuz);
figure(fig)
if(dash)
    plot(delta_vals,nux,'--r');
    hold on
    plot(delta_vals,nuz,'--b');
    xlabel('$\delta$','FontSize',14,'Interpreter','latex');
    grid on
else
    plot(delta_vals,nux,'-r');
    hold on
    plot(delta_vals,nuz,'-b');
     xlabel('$\delta$','FontSize',14,'Interpreter','latex');
    grid on
end