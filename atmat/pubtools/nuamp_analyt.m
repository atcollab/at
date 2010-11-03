function [xdx, xdz, zdz]=nuamp_analyt(ring,sym)
%dnudx=nuamp_analyt(ring,sym)
%looks like it works, still not fully tested.
%B. Nash, 3 Nov. 2010

if nargin==1
   sym=1;
end


[lindat,nu,xi]=atlinopt(ring,0,1:length(ring));


betaxz = cat(1, lindat.beta);
bx = betaxz(:,1);
bz = betaxz(:,2);

mu = cat(1, lindat.mu);
mux = mu(:,1);
muz = mu(:,2);

muxT= pi*nu(1);
muzT= pi*nu(2);

sindex=findcells(ring,'PolynomB');
b3vals=zeros(1,length(sindex));
for q=1:length(sindex)
    L=ring{sindex(q)}.Length;
    b3vals(q)= L*ring{sindex(q)}.PolynomB(3);
end

%now we need to add up dnudj contribution from the sextupoles.
% it involves a double sum over pairs of sextupoles.

xdx = 0;
xdz = 0; 
zdz = 0;

sx = sin(muxT);
s3x = sin(3*muxT);
sxp2z = sin(muxT + 2 * muzT);
sxm2z = sin(muxT - 2*muzT);


for s=1:length(sindex)
    for t=1:length(sindex)
        si = sindex(s);
        ti = sindex(t);
        b = b3vals(s)*b3vals(t);
        c =  sqrt(bx(si)*bx(ti));
        
        cx = cos(abs(mux(si) - mux(ti)) - muxT);
        c3x= cos(3*(abs(mux(si) - mux(ti)) - muxT));
        cxp2z = cos(abs(mux(si) - mux(ti)) + 2*abs((muz(si)-muz(ti))) - (muxT + 2* muzT));
        cxm2z = cos(abs(mux(si) - mux(ti)) - 2*abs((muz(si)-muz(ti))) - (muxT - 2* muzT));
        
       xdx = xdx + b * c ^ 3 * (3*cx/sx + c3x/s3x);
       xdz = xdz + b * c * bz(si) * (-4 * bx(ti) * cx / sx + 2 * bz(ti) * cxp2z / sxp2z - 2 * bz(ti) * cxm2z / sxm2z);
       zdz = zdz + b * c * bz(si) * bz(ti) * (4 * cx / sx + cxp2z / sxp2z + cxm2z / sxm2z);
       
    end
end

xdx = -sym*xdx/(32*pi)/1000;
xdz = -sym*xdz/(32*pi)/1000;
zdz = -sym*zdz/(32*pi)/1000;