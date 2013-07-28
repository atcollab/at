function ring = atloadfielderrs(ring,fielderrstruct)
% ATLOADFIELDERRS will load a field error structure into a ring
% field error structure has the fields: 
% elemindex: element indices in ring to impact
% Nval: reference mpole #, e.g. 2 for Quad, 3 for Sextupole
% nval: multipole index to change
% Bval: Value of relative normal coefficient
% Aval: Value of relative skew coefficient
% radius: reference radius used to compute error

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