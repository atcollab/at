function DiffMat = quantumDiff(elems, varargin)
%quantumDiff    Compute the radiation-diffusion matrix
%
%DIFFMAT=QUANTUMDIFF(RING)
%   RING:       Closed ring AT structure, containing radiative elements and
%               RF cavity. Radiative elements are identified by a
%               PassMethod ending with 'RadPass'.
%
%DIFFMAT=QUANTUMDIFF(RING,RADINDEX)
%   RADINDEX:   Indices of elements where diffusion occurs, typically dipoles
%               and quadrupoles.
%
%DIFFMAT=QUANTUMDIFF(LINE,RADINDEX,ORBIT0)    (Deprecated syntax)
%DIFFMAT=QUANTUMDIFF(...,'orbit',ORBITIN)
%   ORBITIN:    Initial 6-D closed orbit.
%               In this mode, LINE may be a section of a ring.

NumElements=length(elems);

[orb0,varargs]=getoption(varargin, 'orbit', []);
if length(varargs) >= 2	% QUANTUMDIFF(RING,RADINDEX,ORBIT0)
    orb0 = varargs{2};
end
if length(varargs) >= 1 % QUANTUMDIFF(RING,RADINDEX,...)
    radindex=varargs{1};
else                    % QUANTUMDIFF(RING)
    radindex=atgetcells(elems,'PassMethod',@(elem, pass) endsWith(pass,'RadPass'));
end

%[mring, ms, orbit] = findm66(ring,1:NumElements+1);
if isempty(orb0)
    orb=num2cell(findorbit6(elems, 1:NumElements),1)';
else
    orb=num2cell(linepass(elems, orb0, 1:NumElements),1)';
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

