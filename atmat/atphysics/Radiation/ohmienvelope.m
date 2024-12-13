function [envelope, rmsdp, rmsbl, varargout] = ohmienvelope(ring,radindex,refpts) %#ok<INUSD>
%OHMIENVELOPE calculates equilibrium beam envelope in a
% circular accelerator using Ohmi's beam envelope formalism [1].
% [1] K.Ohmi et al. Phys.Rev.E. Vol.49. (1994)
%
% [ENVELOPE, RMSDP, RMSBL] = OHMIENVELOPE(RING,RADELEMINDEX)
% [ENVELOPE, RMSDP, RMSBL] = OHMIENVELOPE(RING,RADELEMINDEX,REFPTS)
%
% RING    - an AT ring.
% RADELEMINDEX - ignored, kept for compatibility
% REFPTS  - reference points along the ring. Default: 1
%
% ENVELOPE is a structure with fields
% Sigma   - [SIGMA(1); SIGMA(2)] - RMS size [m] along
%           the principal axis of a tilted ellips
%           Assuming normal distribution exp(-(Z^2)/(2*SIGMA))
% Tilt    - Tilt angle of the XY ellipse [rad]
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
%
% See also ATENABLE_6D FINDM66

check_6d(ring,true,'strict',0);

NumElements = length(ring);
if nargin<3, refpts=1; end

orb0=findorbit(ring);
[BCUM,Batbeg]=atdiffmat(ring,'orbit',orb0);
[mring, ms, orbit] = findm66(ring,1:NumElements+1,'orbit',orb0);

% ------------------------------------------------------------------------
% Equation for the moment matrix R is
%         R = MRING*R*MRING' + BCUM;
% We rewrite it in the form of Sylvester-Lyapunov equation
% to use MATLAB's SYLVESTER function:
%            AA*R + R*BB = CC
% where
%				AA = inv(MRING)
%				BB = -MRING'
%				CC = inv(MRING)*BCUM
% -----------------------------------------------------------------------
AA = inv(mring);
BB = -mring';
CC = AA*BCUM; %#ok<MINV>

R = sylvester(AA,BB,CC);     % Envelope matrix at the ring entrance

rmsdp = sqrt(R(5,5));   % R.M.S. energy spread
rmsbl = sqrt(R(6,6));   % R.M.S. bunch length

mt=squeeze(num2cell(ms,[1 2]));
[rr,tt,ss]=cellfun(@propag,mt(refpts),Batbeg(refpts),'UniformOutput',false);
envelope=struct('R',rr,'Sigma',ss,'Tilt',tt);

nout=nargout-3;
varargout=cell(1,nout);
if nout>=3, varargout{3}=orbit(:,refpts); end
if nout>=2, varargout{2}=ms(:,:,refpts); end
if nout>=1, varargout{1}=mring; end

    function [r,tilt,sigma]=propag(m,cumb)
        r=m*R*m'+cumb;
        [u,dr] = eig(r([1 3],[1 3]));
        tilt = asin((u(2,1)-u(1,2))/2);
        sigma=sqrt([dr(1,1);dr(2,2)]);
    end

end
