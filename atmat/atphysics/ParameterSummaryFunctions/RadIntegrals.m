function [ I1, I2 ,I3 ,I4 ,I5] = RadIntegrals(ring,Wigidx,Alph,Bta)
%calcuate the contribution to the radiation integrals of a Wiggler. 
%Wigidx is index of insertion devices
%
% A.Mash'al, Iranian Light Source Facility, 2018-07-30   

B0=ring{Wigidx(1),1}.Bmax;
LT=ring{Wigidx(1),1}.Length;
Lw=ring{Wigidx(1),1}.Lw;
kw=2*pi/Lw;
poles=LT/Lw;
alpha(1)=Alph;
beta(1)=Bta;
gamma(1)=(1+alpha(1)^2)/beta(1);
for i=1:poles-1
   beta(i+1)  = beta(i) - 2*Lw*alpha(i) + sqrt(Lw)*gamma(i);
   alpha(i+1) = alpha(i) - Lw*gamma(i);
   gamma(i+1) = (1+alpha(i+1)*alpha(i+1))/beta(i+1);
end

if isempty(ring{Wigidx(1)}.Bx) 
pBy=ring{Wigidx(1),1}.By;
nB=size(pBy,2);
kz=kw*pBy(5,1);
By=@(x) -B0*pBy(2,1)*cos(kz*x+pBy(6,1));
if nB>1
for i=1:nB
    kz=kw*pBy(5,i);
    By0=@(x) -B0*pBy(2,i)*cos(kz*x+pBy(6,i));
    By=By+By0;
    clear('By0')
end
end
end

if isempty(ring{Wigidx(1)}.By)
pBx=ring{Wigidx(1),1}.Bx;
nB=size(pBx,2);
kz=kw*pBx(5,1);
By=@(x) B0*pBx(2,1)*cos(kz*x+pBx(6,1));
if nB>1
for i=1:nB
    kz=kw*pBx(5,i);
    By0=@(x) B0*pBx(2,i)*cos(kz*x+pBx(6,i));
    By=By+By0;
    clear('By0')
end
end
end
helical=0;
if ~isempty(ring{Wigidx(1)}.By) && ~isempty(ring{Wigidx(1)}.Bx)
helical=1;
pBx=ring{Wigidx(1),1}.Bx;
nB=size(pBx,2);
kz=kw*pBx(5,1);
By=@(x) B0*pBx(2,1)*cos(kz*x+pBx(6,1));
if nB>1
for i=1:nB
    kz=kw*pBx(5,i);
    By0=@(x) B0*pBx(2,i)*cos(kz*x+pBx(6,i));
    By=By+By0;
    clear('By0')
end
end
end
     
    
    
    
    
    function dy = Eta(x,y)
    dy = zeros(2,1);    
    dy(1) = y(2);
    dy(2) = -(By(x).*By(x)./100)*y(1)+By(x)/10;
    end
[S,Y] = ode15s(@Eta,[0 LT],[0 0]);
Eta_x=Y(:,1);
Eta_px=Y(:,2);
fI1=Eta_x.*By(S)./10;
I1=trapz(S,fI1);
rho2=@(x) (By(x).^2)./100;
I2=integral(rho2,0,LT);
rho3=@(x) abs(By(x).^3)./1000;
I3=integral(rho3,0,LT);
fI4=Eta_x.*By(S).*(By(S).^2)./1000;
I4=trapz(S,fI4);
alphax=Alph;
betax=Bta;
rhw=10/B0;
gammax=(1+alphax^2)/betax;

I5=(((Lw^4)/(4*(pi^4)*(rhw^5)))*((3/5/pi)+(3/16))*gammax*LT)-...
    (((9*(Lw^3))/(40*(pi^4)*(rhw^5)))*alphax*LT)+...
    ((Lw^2)*betax*LT/(15*(pi^3)*(rhw^5)));

if helical==1
    clear('By')
    pBy=ring{Wigidx(1),1}.By;
nB=size(pBy,2);
kz=kw*pBy(5,1);
By=@(x) -B0*pBy(2,1)*cos(kz*x+pBy(6,1));
if nB>1
for i=1:nB
    kz=kw*pBy(5,i);
    By0=@(x) -B0*pBy(2,i)*cos(kz*x+pBy(6,i));
    By=By+By0;
    clear('By0')
end
end
[S,Y] = ode15s(@Eta,[0 LT],[0 0]);
Eta_x=Y(:,1);
Eta_px=Y(:,2);
fI1=Eta_x.*By(S)./10;
I12=trapz(S,fI1);
rho2=@(x) (By(x).^2)./100;
I22=integral(rho2,0,LT);
rho3=@(x) abs(By(x).^3)./1000;
I32=integral(rho3,0,LT);
fI4=Eta_x.*By(S).*(By(S).^2)./1000;
I42=trapz(S,fI4);
alphax=Alph;
betax=Bta;
gammax=(1+alphax^2)/betax;
H=gammax*(Eta_x.^2)+2*alphax*(Eta_px.^2)+betax*(Eta_x.*Eta_px);
fI5=H.*rho3(S);
I52=trapz(S,fI5);
    I1=I1+I12;
    I2=I2+I22;
    I3=I3+I32;
    I4=I4+I42;
    I5=I5+I52;
end
end
    
