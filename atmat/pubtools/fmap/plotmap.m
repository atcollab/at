function [nuh,nuv]=plotmap(X,Y,nuh,nuv)
%Uses the outputs rom the function "nonlinmap" to calculate the difusion of
%tunes and plot the map

deltnu=log10(sqrt(nuh(:,3)+nuv(:,3)));

matdeltnu=reshape(deltnu,numel(X),numel(Y));
matdeltnu=matdeltnu';
matmeanh=(reshape(nuh(:,4),numel(X),numel(Y)))';
matmeanv=(reshape(nuv(:,4),numel(X),numel(Y)))';
matmean=sqrt(matmeanh.^2+matmeanv.^2);
clf
matmaxh=(reshape(nuh(:,5),numel(X),numel(Y)))';
matmaxv=(reshape(nuv(:,5),numel(X),numel(Y)))';
matmax=sqrt(matmaxh.^2+matmaxv.^2);
%tune difusion versus initial conditions
subplot(2,3,1)
[ax,ay]=meshgrid(X,Y);
hold off
pl=pcolor(ax,ay,matdeltnu);
caxis([-10 0]);
Xlabel('X');
Ylabel('Y');
Title('Tune diffusion');

%tune foot print 
subplot(2,3,2);

mnuh=(reshape(nuh(:,1),numel(X),numel(Y)))';
mnuv=(reshape(nuv(:,1),numel(X),numel(Y)))';

pl=pcolor(mnuh,mnuv,matdeltnu);
set(pl,'FaceColor', 'none','edgecolor','none','marker','.','markeredgecolor','flat','markerfacecolor','flat');
caxis([-10 0]);
Xlabel('Nux');
Ylabel('Nuv');

%Tune shift versus Y amplitude (for X=0 and first X above 0)
subplot(2,3,4);

[Ic]=find((Y>-1E-12)&(Y<1E-12));
[Il]=find((X>-1E-12)&(X<1E-12));


%pl=pcolor(X,Y,matmeanv);
pl=plot(Y'.*Y',mnuv(:,Il)');
%set(pl,'FaceColor', 'none','edgecolor','none','marker','.','markeredgecolor','flat','markerfacecolor','flat');
%caxis([-10 0]);

Xlabel('Y^2 (m*m)');
Ylabel('Nuz');
%Tune shift versus X amplitude ( (for X=0 and first X above 0)
subplot(2,3,3);
%pl=pcolor(X,Y,matmeanh);
pl=plot(X'.*X',mnuh(Ic,:));
Xlabel('X^2 (m*m)');
Ylabel('Nux');
%set(pl,'FaceColor', 'none','edgecolor','none','marker','.','markeredgecolor','flat','markerfacecolor','flat');
%caxis([-10 0]);
subplot(2,3,5);
%pl=pcolor(X,Y,matmeanh);
pl=plot(X'.*X',mnuv(Ic,:));
Xlabel('X^2 (m*m)');
Ylabel('Nuz');
subplot(2,3,6);
%pl=pcolor(X,Y,matmeanh);
pcolor(ax,ay,matmean)
Xlabel('X');
Ylabel('Y');
Title('Mean Orbit')