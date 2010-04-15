function newring=scalesext(ring,sextfam,scale)
newring=ring;
ind=findcells(newring,'FamName',sextfam);
for j=1:length(ind)
    newring{ind(j)}.PolynomB(3)= newring{ind(j)}.PolynomB(3)*scale;
end