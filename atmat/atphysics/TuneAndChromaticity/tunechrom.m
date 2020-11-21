function [tune, varargout] = tunechrom(RING,DP,varargin)
%TUNECHROM computes linear tunes and chromaticities for COUPLED or UNCOUPLED lattice
%
% TUNE = TUNECHROM(RING,DP) - quick calculation of fractional part of the tune 
%    from numerically computed transfer matrix, assuming NO X-Y coupling.
%   
% [TUNE, CHROM] = TUNECHROM(RINGD,DP,'chrom') - optionally computes
%    chromaticity by numerical differentiation from the difference between tune
%    values at momentums DP+0.5*DPStep and DP-0.5*DPStep
%
% [...]=TUNECHROM(..., 'coupling')
% [...]=TUNECHROM(..., 'coupled',boolflag) - Tunes and chromaticities are
%   calculated assuming COUPLED lattice for two transverse eigenmodes.
%
% [...] = TUNECHROM(...,'orbit',ORBITIN - do not search for closed orbit.
%   Instead ORBITIN,a 6x1 vector of initial conditions is used:
%   [x0; px0; y0; py0; DP; 0]. The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit
%   if RING is not a ring or to avoid recomputing the
%   closed orbit if is already known.
%
% Note: TUNECHROM computes tunes and chromaticities from the 4-by-4
%   transfer matrix. The transfer matrix is found in FINDM44 using
%   numerical differentiation. The error of numerical differentiation 
%   is sensitive to the step size. (Reference: Numerical Recipes)
%   Calculation of tunes in TUNECHROM involves one numerical differentiation
%   to find the 4-by-4 transfer matrix.
%   Calculation of chromaticity in TUNECHROM involves TWO!!! numerical differentiations.
%   The error in calculated chromaticity from may be substantial (~ 1e-5).
%   Use the DDP argument to control the step size in chromaticity calculations
%   Another  way to control the step size is NUMDIFPARAMS structure
%   
%   
% See also LINOPT, TWISSRING, TWISSLINE, NUMDIFPARAMS

[COUPLINGFLAG,varargs]=getflag(varargin,'coupling');
[CHROMFLAG,varargs]=getflag(varargs,'chrom');
[orbitin,varargs]=getoption(varargs,'orbit',missing);
[COUPLINGFLAG,varargs]=getoption(varargs,'coupled',COUPLINGFLAG);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');

M44 = findm44(RING,DP,'orbit',orbitin,'XYStep',XYStep);

M =M44(1:2,1:2);
N =M44(3:4,3:4);
if COUPLINGFLAG
    m =M44(1:2,3:4);
    n =M44(3:4,1:2);
    
    % 2-by-2 symplectic matrix
    S = [0 1; -1 0];
    H = m + S*n'*S';
    t = trace(M-N);
    
    g = sqrt(1 + sqrt(t*t/(t*t+4*det(H))))/sqrt(2);
    G = diag([g g]);
    C = -H*sign(t)/(g*sqrt(t*t+4*det(H)));
    A = G*G*M  -  G*(m*S*C'*S' + C*n) + C*N*S*C'*S';
    B = G*G*N  +  G*(S*C'*S'*m + n*C) + S*C'*S'*M*C;
    tune = [closure(A) closure(B)];
else
    tune = [closure(M) closure(N)];
end

if CHROMFLAG
    tune_P = tunechrom(RING,DP+0.5*DPStep,'coupled',COUPLINGFLAG);
    tune_M = tunechrom(RING,DP-0.5*DPStep,'coupled',COUPLINGFLAG);
    varargout{1} = (tune_P - tune_M)/DPStep;
end

    function tune = closure(AB)
        cosmu = (AB(1,1) + AB(2,2))/2;
        diff  = (AB(1,1) - AB(2,2))/2;
        sinmu = sign(AB(1,2))*sqrt(-AB(1,2)*AB(2,1)-diff*diff);
        tune = mod(atan2(sinmu,cosmu)/2/pi,1);
    end        
end
