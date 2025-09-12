function R=VaryQuadFam(R,K1val,fam)
%
% functions returns a new ring with a different value for K1 of fam

indfam=findcells(R,'FamName',fam);
%*ones(size(indfam))
R=setcellstruct(R,'PolynomB',indfam,K1val*ones(size(indfam)),1,2); 

