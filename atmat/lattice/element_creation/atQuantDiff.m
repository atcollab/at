function elem=atQuantDiff(fname,ring)
%   atQuantDiff creates a quantum diffusion element
%   fname is the element name, 
%   ring is the at lattice without radiation
%

[ring2,radindex]=atradon(ring);
dmat=quantumDiff(ring2,radindex);
Lmatp=lmatp(dmat);
%[ld,prm]=atx2(ring);
%diff = prm.DiffCum;
%Lmat=chol(diff);
%Lmatp=Lmat';

elem=atbaselem(fname,'QuantDiffPass','Class','QuantDiff','Lmatp',Lmatp,'Length',0.0);

end
    

function [ lmatp ] = lmatp( dmat )
%lmat does Cholesky decomp of dmat unless diffusion is 0 in
%vertical.  Then do chol on 4x4 hor-long matrix and put 0's
%in vertical

try
    lmat66 = chol(dmat);
catch
    lm=[chol(dmat([1 2 5 6],[1 2 5 6])) zeros(4,2);zeros(2,6)];
    lmat66=lm([1 2 5 6 3 4],[1 2 5 6 3 4]);
end
lmatp=lmat66';
end


