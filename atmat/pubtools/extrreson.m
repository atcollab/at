function res=extrreson(X,Y,nuh,nuv,Hnuh,Hnuv,eps)
deltnu=log10(sqrt(nuh(:,3)+nuv(:,3)));

matdeltnu=reshape(deltnu,numel(X),numel(Y))';

matnuh=(reshape(nuh(:,1),numel(X),numel(Y)))';
matnuv=(reshape(nuv(:,1),numel(X),numel(Y)))';
[ax,ay]=meshgrid(X,Y);

matres=abs(round(matnuh*Hnuh+matnuv*Hnuv)-(matnuh*Hnuh+matnuv*Hnuv));
[I,J]=find(((matres)<eps)|((1-matres)<eps));
size(I)
for i=1:numel(I)
    nx(i)=ax(I(i),J(i));
    ny(i)=ay(I(i),J(i));
    nux(i)=matnuh(I(i),J(i));
    nuz(i)=matnuv(I(i),J(i));
end
%tune difusion versus initial conditions
subplot(2,2,1)
plot(nx,ny,'o');


%tune foot print 
subplot(2,2,2);
plot(nux,nuz,'o');




