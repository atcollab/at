function [BCUM,Batbeg] = atdiffmat(ring,energy,varargin)
%quantumDiff    Compute the radiation-diffusion matrix
%
%[BCUM,BS]=ATDIFFMAT(RING,ENERGY)
%   RING:       Closed ring AT structure, containing radiative elements and
%               RF cavity. Radiative elements are identified by a
%               PassMethod ending with 'RadPass'.
%   ENERGY:     Lattice energy
%
%   BCUM:       Cumulated diffusion matrix
%   BS:         Cumulative diffusion matrix at the beginning of each element
%
%[BCUM,BS]=ATDIFFMAT(...,'orbit',ORBITIN)
%   ORBITIN:    Initial 6-D closed orbit.
%               In this mode, RING may be a section of a ring.

NumElements=length(ring);

[orb0,varargs]=getoption(varargin, 'orbit', []); %#ok<ASGLU>
if isempty(orb0)
    orbit=num2cell(findorbit6(ring, 1:NumElements),1)';
else
    orbit=num2cell(linepass(ring, orb0, 1:NumElements),1)';
end

zr=zeros(6,6);
% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6);
% Batbeg{i} is the cumulative diffusion matrix from
% 0 to the beginning of the i-th element
Batbeg=[zr;cellfun(@cumulb,ring,orbit,'UniformOutput',false)];

    function btx=cumulb(elem,orbit)
        if endsWith(elem.PassMethod,'RadPass')
            if isfield(elem,'Bmax')
                b=FDW(elem,orbit,energy);  % For wigglers
            else
                b=findmpoleraddiffmatrix(elem, orbit, energy);  % For other elements
            end
        else
            b=zr;
        end
        % Calculate 6-by-6 linear transfer matrix in each element
        % near the equilibrium orbit
        m=findelemm66(elem,elem.PassMethod,orbit);
        % Cumulative diffusion matrix of the entire ring
        BCUM = m*BCUM*m' + b;
        btx=BCUM;
    end
end

