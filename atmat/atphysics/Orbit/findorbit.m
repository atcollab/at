function [orbs,orbitin]  = findorbit(ring,varargin)
%FINDORBIT find the closed orbit
%
%Depending on the lattice, FINDORBIT will:
% - use findorbit6 if radiation is ON,
% - use findsyncorbit if radiation is OFF and ct is specified,
% - use findorbit4 otherwise
%
%[ORBIT,FIXEDPOINT]=FINDORBIT(RING,REFPTS)
%
% ORBIT:        6xNrefs, closed orbit at the selected reference points
% FIXEDPOINT:   6x1 closed orbit at the entrance of the lattice
%
%[...]=FINDORBIT(RING,...,'dp',DP) Specify the momentum deviation when
%   radiation is OFF (default: 0)
%
%[...]=FINDORBIT(RING,...,'dct',DCT) Specify the path lengthening when
%   radiation is OFF (default: 0)
%
%[...]=FINDORBIT(RING,...,'orbit',ORBIT) Specify the orbit at the entrance
%   of the ring, if known. FINDORBIT will then transfer it to the reference points
%
%[...]=FINDORBIT(RING,...,'guess',GUESS) The search will start with GUESS as
%   initial conditions
%
%For other keywords, see the underlying functions.
%
% See also FINDORBIT4, FINDSYNCORBIT, FINDORBIT6.

[orbitin,varargs]=getoption(varargin,'orbit',[]);
[refpts,varargs]=getargs(varargs,[],'check',@(arg) isnumeric(arg) || islogical(arg));
if isempty(orbitin)
    if check_radiation(ring)    % Radiation ON: 6D orbit
        dp=getoption(varargs,'dp',NaN);
        ct=getoption(varargs,'dct',NaN);
        if isfinite(dp) || isfinite(ct)
            warning('AT:linopt','In 6D, "dp" and "dct" are ignored');
        end
        orbitin=xorbit_6(ring,varargs{:});
    else                        % Radiation OFF: 4D orbit
        [~,varargs]=getoption(varargs,'DPStep');   % Take it out because findorbit4 does not accept it
        [dp,varargs]=getoption(varargs,'dp',0);
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
    orbs=orbitin;
else
    orbs=linepass(RING,orbitin,refpts,args{:});
end
end
