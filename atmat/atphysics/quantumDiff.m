function [ DiffMat ] = quantumDiff( elems,radindex,orbit )
%quantumDiff gives the diffusion matrix for elems.
%elems may be a ring or section of a ring (or transfer line)
%radindex is the indices for the set of elements where diffusion occurs, typically dipoles
%and quadrupoles.
%orbit is the orbit offset in each element.  If elems represents a closed
%ring, then no orbit input argument is required and the function computes the
%closed orbit.  If elems is not closed, then an array with orbit offset in
%each of the elements in elems is required for the third argument.

NumElements=length(elems);

%[mring, ms, orbit] = findm66(ring,1:NumElements+1);
if (nargin < 3)
    orbit=findorbit6(elems,1:NumElements+1);
end
orb=num2cell(orbit,1)';

zr={zeros(6,6)};
B=zr(ones(NumElements,1));   % B{i} is the diffusion matrix of the i-th element

% calculate Radiation-Diffusion matrix B for elements with radiation
B(radindex)=cellfun(@findmpoleraddiffmatrix,...
    elems(radindex),orb(radindex),'UniformOutput',false);

% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6);
% Batbeg{i} is the cumulative diffusion matrix from
% 0 to the beginning of the i-th element
Batbeg=[zr;cellfun(@cumulb,elems,orb(1:end-1),B,'UniformOutput',false)]; %#ok<NASGU>

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

