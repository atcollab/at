function [ I1, I2 ,I3 ,I4 ,I5] = RadIntegrals(ring,Wigidx,Alph,Bta,Dis0)
%calcuate the contribution to the radiation integrals of a Wiggler. 
%Wigidx is index of insertion devices
%
%--- Modification Log --------------------------------
% A.Mash'al, Iranian Light Source Facility, 2020-09-03
% A.Mash'al, Iranian Light Source Facility, 2018-07-30
B0=ring{Wigidx(1),1}.Bmax;
LT=ring{Wigidx(1),1}.Length;
Lw=ring{Wigidx(1),1}.Lw;
E0=ring{Wigidx(1),1}.Energy/1e9;  %Energy in [GeV]
kw=2*pi/Lw;
%------------ On-axis Magnetic Field B(0,0,z)----------
% Y. Wu,"SYMPLECTIC MODELS FOR GENERAL INSERTION DEVICES", Eq (8) 
% M. Borland, GWigB function @ gwigR.c

% Horizontal wiggler , Bx=[] 
if isempty(ring{Wigidx(1)}.Bx) 
    
WMode='Horz';    
    
pBy=ring{Wigidx(1),1}.By;
nHarm=size(pBy,2);
By=@(x) 0;
for i=1:nHarm
    kz=kw*pBy(5,i);     %k_zn= n*k_w  , n=pBy(5,i)
    tz=pBy(6,i);
    C_mn=pBy(2,i);
    By0=@(x) -B0*C_mn*cos(kz*x+tz);
    By=@(x) By(x)+By0(x);
end
BT=By;
end


% Vertical wiggler , By=[] 
if isempty(ring{Wigidx(1)}.By) 
    
WMode='Vert';
    
pBx=ring{Wigidx(1),1}.Bx;
nHarm=size(pBx,2);
Bx=@(x) 0;
for i=1:nHarm
    kz=kw*pBx(5,i);     %k_zn= n*k_w  , n=pBx(5,i)
    tz=pBx(6,i);
    C_mn=pBx(2,i);
    Bx0=@(x) B0*C_mn*cos(kz*x+tz);
    Bx=@(x) Bx(x)+Bx0(x);
end
BT=Bx;
end

% Elliptical Polarized Wiggler

if ~isempty(ring{Wigidx(1)}.By) && ~isempty(ring{Wigidx(1)}.Bx)

WMode='Ellip';

pBy=ring{Wigidx(1),1}.By;
pBx=ring{Wigidx(1),1}.Bx;
nHarmH=size(pBy,2);
nHarmV=size(pBx,2);
By=@(x) 0;
Bx=@(x) 0;
for i=1:nHarmH
    kz=kw*pBy(5,i);     
    tz=pBy(6,i);
    C_mn=pBy(2,i);
    By0=@(x) -B0*C_mn*cos(kz*x+tz);
    By=@(x) By(x)+By0(x);
end
for i=1:nHarmV
    kz=kw*pBx(5,i);     
    tz=pBx(6,i);
    C_mn=pBx(2,i);
    Bx0=@(x) B0*C_mn*cos(kz*x+tz);
    Bx=@(x) Bx(x)+Bx0(x);
end
BT=@(x) sqrt((Bx(x).^2)+(By(x).^2));
end   
%------------- Dispersion Through Wiggler ----------------
ring{Wigidx(1),1}.Nstep=100;  %increase integration precision
ring{Wigidx(1),1}.Nmeth=4;
WIGELEM=ring(Wigidx(1));
d=1e-3;
Oplus =[ d*Dis0(1);  d*Dis0(2);  d*Dis0(3);  d*Dis0(4);  d;0];
Ominus=[-d*Dis0(1); -d*Dis0(2); -d*Dis0(3); -d*Dis0(4); -d;0];
linepass(WIGELEM,Oplus,1:2);
ORBp = importdata('Diswig.DAT');
sp=ORBp(:,5)';
ORBp=ORBp(:,1:4);
linepass(WIGELEM,Ominus,1:2);
ORBm = importdata('Diswig.DAT');
ORBm=ORBm(:,1:4);
DisW=(ORBp-ORBm)/d/2;
Dsx=DisW(:,1)';
Dspx=DisW(:,2)';
Dsy=DisW(:,3)';
Dspy=DisW(:,4)';
%------------- Radiation Integrals------------------------
rho2=@(x) (BT(x).^2)*((0.2998/E0)^2);
I2=integral(rho2,0,LT);
rho3=@(x) abs(BT(x).^3)*((0.2998/E0)^3);
I3=integral(rho3,0,LT);

switch WMode
    case 'Horz'
        RHOy=@(x) By(x)*0.2998/E0;
        fI1=Dsx.*RHOy(sp);
        I1=trapz(sp,fI1);
        fI4=Dsx.*By(sp).*(By(sp).^2)*((0.2998/E0)^3);
        I4=trapz(sp,fI4);
        alphax=Alph(1)*ones(1,length(sp));
        betax=Bta(1)*ones(1,length(sp));
        gammax=(1+alphax.^2)./betax;
        H=gammax.*(Dsx.^2)+2*alphax.*(Dspx.^2)+betax.*(Dsx.*Dspx);
        fI5=H.*abs(By(sp).^3)*((0.2998/E0)^3);
        I5=trapz(sp,fI5);
    case 'Vert'
        RHOx=@(x) Bx(x)*0.2998/E0;
        fI1=Dsy.*RHOx(sp);
        I1=trapz(sp,fI1);
        fI4=Dsy.*Bx(sp).*(Bx(sp).^2)*((0.2998/E0)^3);
        I4=trapz(sp,fI4);
        alphax=Alph(2)*ones(1,length(sp));
        betax=Bta(2)*ones(1,length(sp));
        gammax=(1+alphax.^2)./betax;
        H=gammax.*(Dsy.^2)+2*alphax.*(Dspy.^2)+betax.*(Dsy.*Dspy);
        fI5=H.*abs(Bx(sp).^3)*((0.2998/E0)^3);
        I5=trapz(sp,fI5);
    case 'Ellip'
        RHOy=@(x) By(x)*0.2998/E0;
        fI1=Dsx.*RHOy(sp);
        I1y=trapz(sp,fI1);
        fI4=Dsx.*By(sp).*(By(sp).^2)*((0.2998/E0)^3);
        I4y=trapz(sp,fI4);
        alphax=Alph(1)*ones(1,length(sp));
        betax=Bta(1)*ones(1,length(sp));
        gammax=(1+alphax.^2)./betax;
        H=gammax.*(Dsx.^2)+2*alphax.*(Dspx.^2)+betax.*(Dsx.*Dspx);
        fI5=H.*abs(By(sp).^3)*((0.2998/E0)^3);
        I5y=trapz(sp,fI5);
        
        RHOx=@(x) Bx(x)*0.2998/E0;
        fI1=Dsy.*RHOx(sp);
        I1x=trapz(sp,fI1);
        fI4=Dsy.*Bx(sp).*(Bx(sp).^2)*((0.2998/E0)^3);
        I4x=trapz(sp,fI4);
        alphax=Alph(2)*ones(1,length(sp));
        betax=Bta(2)*ones(1,length(sp));
        gammax=(1+alphax.^2)./betax;
        H=gammax.*(Dsy.^2)+2*alphax.*(Dspy.^2)+betax.*(Dsy.*Dspy);
        fI5=H.*abs(Bx(sp).^3)*((0.2998/E0)^3);
        I5x=trapz(sp,fI5);
        
        I1=I1y+I1x;
        I4=I4y+I4x;
        I5=I5y+I5x;
end
end