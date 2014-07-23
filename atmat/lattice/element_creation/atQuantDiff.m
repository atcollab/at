function elem=atQuantDiff(fname,ring)
%   atQuantDiff creates a quantum diffusion element
%   fname is the element name, 
%   ring is the at lattice without radiation
%

dmat=quantumDiff(ring);
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
    dmat44=dmat([1 2 5 6],[1 2 5 6]);
    lmat44=chol(dmat44);
    lmat46=[lmat44([1 2],:);zeros(2,4);lmat44([3 4],:)];
    lmat66=[lmat46(:,[1 2]),zeros(6,2),lmat46(:,[3,4])];
end
lmatp=lmat66';
end


function [ DiffMat ] = quantumDiff( ring )
%quantumDiff gives a random kick to simulate the
%global (one turn) effect of quantum diffusion on the electron

%refpts=1:length(ring);

%[ring2,radindex,cavindex,energy]=atradon(ring);
%[envelope,espread,blength,diff,m,T]=ohmienvelope2(ring2,radindex,refpts);

%Lmat=chol(diff);

[ring2,radindex,cavindex,energy]=atradon(ring);
%ring2=atsetcavity(ring2,5e6,1,992);
%ring=ring2;
%ring2=atsetcavity(ring2,RFV,1,harm)

%Compute the diffusion matrix
NumElements=length(ring);

[mring, ms, orbit] = findm66(ring,1:NumElements+1);
mt=squeeze(num2cell(ms,[1 2]));
orb=num2cell(orbit,1)';

zr={zeros(6,6)};
B=zr(ones(NumElements,1));   % B{i} is the diffusion matrix of the i-th element

% calculate Radiation-Diffusion matrix B for elements with radiation
B(radindex)=cellfun(@findmpoleraddiffmatrix,...
    ring(radindex),orb(radindex),'UniformOutput',false);

% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6);
% Batbeg{i} is the cumulative diffusion matrix from
% 0 to the beginning of the i-th element
Batbeg=[zr;cellfun(@cumulb,ring,orb(1:end-1),B,'UniformOutput',false)];

DiffCum = BCUM;

DiffMat=(DiffCum+DiffCum')/2;

%Lmat=chol((DiffCum+DiffCum')/2);

 function btx=cumulb(elem,orbit,b)
        % Calculate 6-by-6 linear transfer matrix in each element
        % near the equilibrium orbit
        m=findelemm66(elem,elem.PassMethod,orbit);
        % Cumulative diffusion matrix of the entire ring
        BCUM = m*BCUM*m' + b;
        btx=BCUM;
 end
end

