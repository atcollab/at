function [ DiffMat ] = quantumDiff( ring,radindex )
%quantumDiff gives a random kick to simulate the
%global (one turn) effect of quantum diffusion on the electron

NumElements=length(ring);

%[mring, ms, orbit] = findm66(ring,1:NumElements+1);
orbit=findorbit6(ring,1:NumElements+1);
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
Batbeg=[zr;cellfun(@cumulb,ring,orb(1:end-1),B,'UniformOutput',false)]; %#ok<NASGU>

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

