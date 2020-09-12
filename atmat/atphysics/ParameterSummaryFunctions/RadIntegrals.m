function [ I1, I2 ,I3 ,I4 ,I5] = RadIntegrals(ring,Wigidx,Alph,Bta,Dis0)
%calcuate the contribution to the radiation integrals of a Wiggler. 
%Wigidx is index of insertion devices
%
%--- Modification Log --------------------------------
% A.Mash'al, Iranian Light Source Facility, 2020-09-03
% A.Mash'al, Iranian Light Source Facility, 2018-07-30

global GLOBVAL
B0=ring{Wigidx(1),1}.Bmax;
LT=ring{Wigidx(1),1}.Length;
Lw=ring{Wigidx(1),1}.Lw;
kw=2*pi/Lw;
poles=LT/Lw;
if isempty(ring{Wigidx(1)}.Bx) 
pBy=ring{Wigidx(1),1}.By;
nB=size(pBy,2);
kz=kw*pBy(5,1);
By=@(x) -B0*pBy(2,1)*sin(kz*x+pBy(6,1));
if nB>1
for i=1:nB
    kz=kw*pBy(5,i);
    By0=@(x) -B0*pBy(2,i)*sin(kz*x+pBy(6,i));
    By=By+By0;
    clear('By0')
end
end
end

if isempty(ring{Wigidx(1)}.By)
pBx=ring{Wigidx(1),1}.Bx;
nB=size(pBx,2);
kz=kw*pBx(5,1);
By=@(x) B0*pBx(2,1)*sin(kz*x+pBx(6,1));
if nB>1
for i=1:nB
    kz=kw*pBx(5,i);
    By0=@(x) B0*pBx(2,i)*sin(kz*x+pBx(6,i));
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
By=@(x) B0*pBx(2,1)*sin(kz*x+pBx(6,1));
if nB>1
for i=1:nB
    kz=kw*pBx(5,i);
    By0=@(x) B0*pBx(2,i)*sin(kz*x+pBx(6,i));
    By=By+By0;
    clear('By0')
end
end
end
E0=GLOBVAL.E0/1e9;
rho2=@(x) (By(x).^2)*((0.2998/E0)^2);
I2=integral(rho2,0,LT);
rho3=@(x) abs(By(x).^3)*((0.2998/E0)^3);
I3=integral(rho3,0,LT);
RHO=@(x) By(x)*0.2998/E0;
close all
[Ds,Dsp,sp]= Dispwig (RHO,Lw,Dis0);
Ds=Ds';
Dsp=Dsp';
fI1=Ds.*By(sp)*0.2998/E0;
I1=poles*trapz(sp,fI1);
disp(I1)
fI4=Ds.*By(sp).*(By(sp).^2)*((0.2998/E0)^3);
I4=poles*trapz(sp,fI4);
alphax=Alph(1)*ones(1,length(sp));
betax=Bta(1)*ones(1,length(sp));
gammax=(1+alphax.^2)./betax;
H=gammax.*(Ds.^2)+2*alphax.*(Dsp.^2)+betax.*(Ds.*Dsp);
fI5=H.*rho3(sp);
I5=poles*trapz(sp,fI5);


if helical==1
    clear('By')
    pBy=ring{Wigidx(1),1}.By;
nB=size(pBy,2);
kz=kw*pBy(5,1);
By=@(x) -B0*pBy(2,1)*sin(kz*x+pBy(6,1));
if nB>1
for i=1:nB
    kz=kw*pBy(5,i);
    By0=@(x) -B0*pBy(2,i)*sin(kz*x+pBy(6,i));
    By=By+By0;
    clear('By0')
end
end
rho2=@(x) (By(x).^2)*((0.2998/E0)^2);
I22=integral(rho2,0,LT);
rho3=@(x) abs(By(x).^3)*((0.2998/E0)^3);
I32=integral(rho3,0,LT);
RHO=@(x) By(x)*0.2998/E0;
[Ds,Dsp,sp]= Dispwig (RHO,Lw,0);
Ds=Ds';
Dsp=Dsp';
fI1=Ds.*By(sp)*0.2998/E0;
I12=poles*trapz(sp,fI1);
disp(I12)
fI4=Ds.*By(sp).*(By(sp).^2)*((0.2998/E0)^3);
I42=poles*trapz(sp,fI4);
alphax=Alph(2)*ones(1,length(sp));
betax=Bta(2)*ones(1,length(sp));
gammax=(1+alphax.^2)./betax;
H=gammax.*(Ds.^2)+2*alphax.*(Dsp.^2)+betax.*(Ds.*Dsp);
fI5=H.*rho3(sp);
I52=poles*trapz(sp,fI5);
    I1=I1+I12;
    I2=I2+I22;
    I3=I3+I32;
    I4=I4+I42;
    I5=I5+I52;
end
end