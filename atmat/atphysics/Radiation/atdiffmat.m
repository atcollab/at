function [BCUM,Batbeg] = atdiffmat(ring,varargin)
%quantumDiff    Compute the radiation-diffusion matrix
%
%[BCUM,BS]=ATDIFFMAT(RING)
%   RING:       Closed ring AT structure, containing radiative elements and
%               RF cavity. Radiative elements are identified by a
%               PassMethod ending with 'RadPass'.
%
%   BCUM:       Cumulated diffusion matrix
%   BS:         Cumulative diffusion matrix at the beginning of each element
%
%[BCUM,BS]=ATDIFFMAT(RING,'orbit',ORBITIN)
%   ORBITIN:    Initial 6-D closed orbit.
%               In this mode, RING may be a section of a ring.

newmethods = {'BndMPoleSymplectic4RadPass', ...
              'StrMPoleSymplectic4RadPass', ...
              'ExactMultipoleRadPass',...
              'GWigSymplecticRadPass',...
              'EnergyLossRadPass'};

NumElements=length(ring);

[orb0,varargs]=getoption(varargin, 'orbit', []); %#ok<ASGLU>

test_mode = getoption('test_mode');
[energy,cell_length,particle]=atGetRingProperties(ring,'energy','cell_length','particle');

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
        passm = elem.PassMethod;
        if ~endsWith(passm, 'RadPass')
            bdiff = zr;
        elseif ~test_mode
            % New method
            bdiff = diffusion_matrix(substitute(elem),orbit,energy,particle,cell_length,0.0);
        else
            % Old method
            if isfield(elem, 'Bmax')
                bdiff = FDW(elem, orbit, energy);
            else
                bdiff = findmpoleraddiffmatrix(elem,orbit,energy);
            end
        end
        % Calculate 6-by-6 linear transfer matrix in each element
        % near the equilibrium orbit
        m=findelemm66(elem,passm,orbit);
        % Cumulative diffusion matrix of the entire ring
        BCUM = m*BCUM*m' + bdiff;
        btx=BCUM;
    end
    
    function elem=substitute(elem)
        if ~any(strcmp(elem.PassMethod, newmethods))
            elem.PassMethod = 'BndMPoleSymplectic4RadPass';
        end
    end
end

