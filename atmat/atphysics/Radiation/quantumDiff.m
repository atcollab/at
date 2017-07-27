function DiffMat = quantumDiff(elems,radindex,orb0)
%quantumDiff    Compute the radiation-diffusion matrix
%
%DIFFMAT=QUANTUMDIFF(RING,RADINDEX)
%
%   RING:       Closed ring AT structure, containing radiative elements and
%               RF cavity
%   RADINDEX:   Indices of elements where diffusion occurs, typically dipoles
%               and quadrupoles.
%
%DIFFMAT=QUANTUMDIFF(ELEMS,RADINDEX,ORBIT0)
%
%   ELEMS:      AT structure
%   RADINDEX:   Indices of elements where diffusion occurs, typically dipoles
%               and quadrupoles.
%   ORBIT0:     Initial 6-D closed orbit
%
% In this mode, ELEMS may be a section of a ring

NumElements=length(elems);

%[mring, ms, orbit] = findm66(ring,1:NumElements+1);
if (nargin >= 3)
    orb=num2cell(linepass(elems, orb0, 1:NumElements),1)';
else
    orb=num2cell(findorbit6(elems, 1:NumElements),1)';
end

zr={zeros(6,6)};
B=zr(ones(NumElements,1));   % B{i} is the diffusion matrix of the i-th element

% calculate Radiation-Diffusion matrix B for elements with radiation
B(radindex)=cellfun(@findmpoleraddiffmatrix,...
    elems(radindex),orb(radindex),'UniformOutput',false);

% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6);
% Batbeg{i} is the cumulative diffusion matrix from
% 0 to the beginning of the i-th element
Batbeg=[zr;cellfun(@cumulb,elems,orb,B,'UniformOutput',false)]; %#ok<NASGU>

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

