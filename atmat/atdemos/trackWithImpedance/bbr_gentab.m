function [ s,bbrx,bbry,bbrz ] = bbr_gentab(xr,npoints,freqx,freqy,freqz,qx,qy,qz,rx,ry,rz)
%xr is extent of wake function in m
%npoints is number of points in table
%freq in GHz
%Rshunt in MOhm

%Based on equations 2.84 and 2.88 in Chao's 'Physics of Collective Instabilities'

clight = 2.99792458e8;
if mod(npoints,2) == 0
    npoints=npoints+1;
end

s = linspace(-xr,xr,npoints);
    
omegarx = 2.0*pi*freqx*1.0e9;
alphax = omegarx/(2.0*qx);
omegabarx = sqrt(omegarx*omegarx-alphax*alphax);

omegary = 2.0*pi*freqy*1.0e9;
alphay = omegary/(2.0*qy);
omegabary = sqrt(omegary*omegary-alphay*alphay);

omegarz = 2.0*pi*freqz*1.0e9;
alphaz = omegarz/(2.0*qz);
omegabarz = sqrt(omegarz*omegarz-alphaz*alphaz);

bbrx=zeros(length(s),1);
bbry=zeros(length(s),1);
bbrz=zeros(length(s),1);

for i = 1:length(s)
    %chao 2.88
    %chao 2.84
    if s(i)==0.0
        bbrz(i)=rz*1e6*alphaz;
        bbrx(i)=0.0;
        bbry(i)=0.0;    
    elseif s(i)<0.0
        bbrx(i)=0.0;
        bbry(i)=0.0;
        bbrz(i)=0.0;
    else
        bbrx(i) = rx*1e6*omegarx*omegarx/qx/omegabarx*exp(alphax*(-s(i))/clight)*sin(omegabarx*(-s(i))/clight);
        bbry(i) = ry*1e6*omegary*omegary/qy/omegabary*exp(alphay*(-s(i))/clight)*sin(omegabary*(-s(i))/clight);
        bbrz(i) = rz*1e6*2.0*alphaz*exp(alphaz*(-s(i))/clight)*(cos(omegabarz*(-s(i))/clight) + ...
               alphaz/omegabarz*sin(omegabarz*(-s(i))/clight)); 
    end
end

end