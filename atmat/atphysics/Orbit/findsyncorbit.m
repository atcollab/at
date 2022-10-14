function [orb5,orbitin] = findsyncorbit(ring,varargin)
%FINDSYNCORBIT finds closed orbit, synchronous with the RF cavity
% and momentum deviation dP (first 5 components of the phase space vector)
% by numerically solving  for a fixed point
% of the one turn map M calculated with LINEPASS
%
%       (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
%
%    under constraints CT2 - CT1 =  dCT = C(1/Frev - 1/Frev0) and dP2 = dP1 , where
%    Frev0 = Frf0/HarmNumber is the design revolution frequency
%    Frev  = (Frf0 + dFrf)/HarmNumber is the imposed revolution frequency
%
% IMPORTANT!!!  FINDSYNCORBIT imposes a constraint (CT2 - CT1) and
%    dP2 = dP1 but no constraint on the value of dP1, dP2
%    The algorithm assumes time-independent fixed-momentum ring
%    to reduce the dimensionality of the problem.
%    To impose this artificial constraint in FINDSYNCORBIT
%    PassMethod used for any element SHOULD NOT
%    1. change the longitudinal momentum dP (cavities , magnets with radiation)
%    2. have any time dependence (localized impedance, fast kickers etc).
%
%
% FINDSYNCORBIT(RING) is 5x1 vector - fixed point at the
%		entrance of the 1-st element of the RING (x,px,y,py,dP)
%
% FINDSYNCORBIT(RING,REFPTS) is 5xLength(REFPTS)
%   array of column vectors - fixed points (x,px,y,py,dP)
%   at the entrance of each element indexed by the REFPTS array.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(RING)+1.
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
%FINDORBIT4(RING,DCT,REFPTS,...) Obsolete syntax
%FINDORBIT4(RING,...,'dct',DCT)  Specifies the path lengthening
%
%FINDORBIT4(RING,...,'dp',DP)   Specifies the off-momentum
%
%FINDORBIT4(RING,...,'df',DF) Specifies RF frequency deviation
%
% FINDSYNCORBIT(RING,DCT,REFPTS,GUESS)
% FINDSYNCORBIT(...,'guess',GUESS)     The search for the fixed point
%   starts from initial condition GUESS. Otherwise the search starts from
%   [0; 0; 0; 0; 0; 0]. GUESS must be a 6x1 vector.
%
% FINDSYNCORBIT(...,'orbit',ORBIT)     Specify the orbit at the entrance
%   of the ring, if known. FINDSYNCORBIT will then transfer it to the
%   reference points. ORBIT must be a 6x1 vector.
%
% [ORBIT, FIXEDPOINT] = FINDSYNCORBIT( ... )
%	The optional second return parameter is a 6x1 vector:
%   closed orbit at the entrance of the RING.
%
% See also FINDORBIT4, FINDORBIT6.
%
if ~iscell(ring)
    error('First argument must be a cell array');
end
[orbitin,varargs]=getoption(varargin,'orbit',[]);
[dct,varargs]=getdparg(varargs,0.0,'key','dct');
[dp,varargs]=getoption(varargs,'dp',NaN);
[df,varargs]=getoption(varargs,'df',NaN);
[refpts,varargs]=getargs(varargs,[],'check',@(arg) isnumeric(arg) || islogical(arg));
if isempty(orbitin)
    if isfinite(df)
        [cell_l,cell_frev,cell_h]=atGetRingProperties(ring,'cell_length',...
            'cell_harmnumber','cell_revolution_frequency');
        dct=-cell_l*df/(cell_frev*cell_h+df);
        orbitin=xorbit_ct(ring,dct,varargs{:});
    elseif isfinite(dp)
        orbitin=xorbit_dp(ring,dp,varargs{:});
    else
        orbitin=xorbit_ct(ring,dct,varargs{:});
    end
    args={'KeepLattice'};
else
    args={};
end

if islogical(refpts)
    refpts=find(refpts);
end
if isempty(refpts)
    % return only the fixed point at the entrance of RING{1}
    orb5=orbitin(1:5,1);
else
    orb6 = linepass(ring,orbitin,refpts,args{:});
    orb5 = orb6(1:5,:);
end
end
