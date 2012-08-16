function ring = atloadfielderrs(ring,fielderrstruct)
% ATLOADFIELDERRS will load a field error structure into a ring

elemindex=fielderrstruct.elemindex;
newindex=fielderrstruct.nval;
refindex=fielderrstruct.Nval;
Bval=fielderrstruct.Bval;
Aval=fielderrstruct.Aval;
radius=fielderrstruct.radius;

for j=1:length(elemindex);
    elem=ring{elemindex(j)};
    polyB=elem.PolynomB;
    polyA=elem.PolynomA;
    
    
    polyB=ataddmpolecomppoly(polyB,refindex(j),newindex(j),Bval(j),radius);
    polyA=ataddmpolecomppoly(polyA,refindex(j),newindex(j),Aval(j),radius,polyB(refindex(j)));
    
    
    elem.PolynomB=polyB;
    elem.PolynomA=polyA;
    elem.MaxOrder=length(polyB);
    ring{elemindex(j)}=elem;
end