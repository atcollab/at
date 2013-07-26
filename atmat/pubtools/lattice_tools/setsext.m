function newring=setsext(ring,sextfam,val)
%newring=setsext(ring,fam,val);
%note: use integrated strengths.  If it is a thick sext, then the length
%will be divided out to put the correct value in the polynomB
newring=ring;
ind=findcells(newring,'FamName',sextfam);
   
for j=1:length(ind)
    val0=val;
    if strcmp(ring{ind(j)}.BetaCode,'SX'), val0=val0/(ring{ind(j)}.Length); end;
    newring{ind(j)}.PolynomB(3)= val0;
end