function [orb4,orbitin] = findorbit4(ring,varargin)
%FINDORBIT4 finds closed orbit in the 4-d transverse phase
% space by numerically solving  for a fixed point of the one turn
% map M calculated with LINEPASS
%
%         (X, PX, Y, PY, dP2, CT2 ) = M (X, PX, Y, PY, dP1, CT1)
%
%    under the CONSTANT MOMENTUM constraint, dP2 = dP1 = dP and
%    there is NO constraint on the 6-th coordinate CT
%
% IMPORTANT!!! FINDORBIT4 imposes a constraint on dP and relaxes
%    the constraint on the revolution frequency. A physical storage
%    ring does exactly the opposite: the momentum deviation of a
%    particle on the closed orbit settles at the value
%    such that the revolution is synchronous with the RF cavity
%
%                 HarmNumber*Frev = Frf
%
%    To impose this artificial constraint in FINDORBIT4
%    PassMethod used for any elemen SHOULD NOT
%    1. change the longitudinal momentum dP (cavities , magnets with radiation)
%    2. have any time dependence (localized impedance, fast kickers etc)
%
%FINDORBIT4(RING) is 4x1 vector - fixed point at the
%    entrance of the 1-st element of the RING (x,px,y,py)
%
%FINDORBIT4(RING,REFPTS) is 4xLength(REFPTS)
%   array of column vectors - fixed points (x,px,y,py)
%   at the entrance of each element indexed by the REFPTS array.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(RING)+1.
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
%FINDORBIT4(RING,DP,REFPTS,...) Obsolete syntax
%FINDORBIT4(RING,...,'dp',DP)   Specifies the off-momentum
%
%FINDORBIT4(RING,...,'dct',DCT) Specifies the path lengthening
%
%FINDORBIT4(RING,...,'df',DF) Specifies RF frequency deviation
%
%FINDORBIT4(RING,dP,REFPTS,GUESS)
%FINDORBIT4(...,'guess',GUESS)     The search for the fixed point
%   starts from initial condition GUESS. Otherwise the search starts from
%   [0; 0; 0; 0; 0; 0]. GUESS must be a 6x1 vector.
%
%FINDORBIT4(...,'orbit',ORBIT)     Specify the orbit at the entrance
%   of the ring, if known. FINDORBIT4 will then transfer it to the
%   reference points. ORBIT must be a 6x1 vector.
%
%[ORBIT, FIXEDPOINT] = FINDORBIT4( ... )
%	The optional second return parameter is a 6x1 vector:
%   closed orbit at the entrance of the RING.
%
% See also FINDSYNCORBIT, FINDORBIT6.

if ~iscell(ring)
    error('First argument must be a cell array');
end
[orbitin,varargs]=getoption(varargin,'orbit',[]);
[dp,varargs]=getdparg(varargs,0.0);
[dct,varargs]=getoption(varargs,'dct',NaN);
[df,varargs]=getoption(varargs,'df',NaN);
[refpts,varargs]=getargs(varargs,[],'check',@(arg) isnumeric(arg) || islogical(arg));
[~,varargs]=getoption(varargs,'is_6d',[]); % Consume the is_6d option
if isempty(orbitin)
    if isfinite(df)
        [cell_l,cell_frev,cell_h]=atGetRingProperties(ring,'cell_length',...
            'cell_harmnumber','cell_revolution_frequency');
        dct=-cell_l*df/(cell_frev*cell_h+df);
        orbitin=xorbit_ct(ring,dct,varargs{:});
    elseif isfinite(dct)
        orbitin=xorbit_ct(ring,dct,varargs{:});
    else
        orbitin=xorbit_dp(ring,dp,varargs{:});
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
    orb4=orbitin(1:4,1);
else
    orb6 = linepass(ring,orbitin,refpts,args{:});
    orb4 = orb6(1:4,:);
end
end