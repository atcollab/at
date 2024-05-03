function [envelope, rmsdp, rmsbl, varargout] = ohmienvelope(ring,radindex,refpts)
%OHMIENVELOPE calculates equilibrium beam envelope in a
% circular accelerator using Ohmi's beam envelope formalism [1]
% [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
%
% [ENVELOPE, RMSDP, RMSBL] = ONMIENVELOPE(RING,RADELEMINDEX,REFPTS)
%
% ENVELOPE is a structure with fields
% Sigma   - [SIGMA(1); SIGMA(2)] - RMS size [m] along
%           the principal axis of a tilted ellips
%           Assuming normal distribution exp(-(Z^2)/(2*SIGMA))
% Tilt    - Tilt angle of the XY ellips [rad]
%           Positive Tilt corresponds to Corkscrew (right)
%           rotatiom of XY plane around s-axis
% R       - 6-by-6 equilibrium envelope matrix R
%
% RMSDP   - RMS momentum spread
% RMSBL   - RMS bunch length[m]
%
% [ENVELOPE, RMSDP, RMSBL, M66, T, ORBIT] = OHMIENVELOPE(...)
%   Returns in addition the 6x6 transfer matrices and the closed orbit
%   from FINDM66

NumElements = length(ring);
if nargin<3, refpts=1; end

Wig=atgetcells(ring,'Bmax');
Wigidx=find(Wig(:)==1);

radindex(Wigidx) = false; % Erase wigglers from the radiative element list.
% Diffusion matrix to be computed with separate FDW function.

[mring, ms, orbit] = findm66(ring,1:NumElements+1);
mt=squeeze(num2cell(ms,[1 2]));
orb=num2cell(orbit,1)';

zr={zeros(6,6)};
B=zr(ones(NumElements,1));   % B{i} is the diffusion matrix of the i-th element

% calculate Radiation-Diffusion matrix B for elements with radiation
B(radindex)=cellfun(@findmpoleraddiffmatrix,...
    ring(radindex),orb(radindex),'UniformOutput',false);
B(Wig)=cellfun(@FDW,...
    ring(Wig),orb(Wig),'UniformOutput',false);

% Calculate cumulative Radiation-Diffusion matrix for the ring
BCUM = zeros(6,6);
% Batbeg{i} is the cumulative diffusion matrix from
% 0 to the beginning of the i-th element
Batbeg=[zr;cellfun(@cumulb,ring,orb(1:end-1),B,'UniformOutput',false)];

% ------------------------------------------------------------------------
% Equation for the moment matrix R is
%         R = MRING*R*MRING' + BCUM;
% We rewrite it in the form of Sylvester-Lyapunov equation
% to use MATLAB's SYLVERTER function:
%            AA*R + R*BB = CC
% where
%				AA = inv(MRING)
%				BB = -MRING'
%				CC = inv(MRING)*BCUM
% -----------------------------------------------------------------------
AA = inv(mring);
BB = -mring';
CC = AA*BCUM;

R = sylvester(AA,BB,CC);     % Envelope matrix at the ring entrance

rmsdp = sqrt(R(5,5));   % R.M.S. energy spread
rmsbl = sqrt(R(6,6));   % R.M.S. bunch length

[rr,tt,ss]=cellfun(@propag,mt(refpts),Batbeg(refpts),'UniformOutput',false);
envelope=struct('R',rr,'Sigma',ss,'Tilt',tt);

nout=nargout-3;
varargout=cell(1,nout);
if nout>=3, varargout{3}=orbit(:,refpts); end
if nout>=2, varargout{2}=ms(:,:,refpts); end
if nout>=1, varargout{1}=mring; end

    function btx=cumulb(elem,orbit,b)
        % Calculate 6-by-6 linear transfer matrix in each element
        % near the equilibrium orbit
        m=findelemm66(elem,elem.PassMethod,orbit);
        % Cumulative diffusion matrix of the entire ring
        BCUM = m*BCUM*m' + b;
        btx=BCUM;
    end

    function [r,tilt,sigma]=propag(m,cumb)
        r=m*R*m'+cumb;
        [u,dr] = eig(r([1 3],[1 3]));
        tilt = asin((u(2,1)-u(1,2))/2);
        sigma=sqrt([dr(1,1);dr(2,2)]);
    end

end
