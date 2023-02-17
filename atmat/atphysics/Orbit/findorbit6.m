function [orb6,orbitin] = findorbit6(ring,varargin)
%FINDORBIT6 finds closed orbit in the full 6-d phase space
% by numerically solving  for a fixed point of the one turn
% map M calculated with RINGPASS
%
% (X, PX, Y, PY, DP, CT2 ) = M (X, PX, Y, PY, DP, CT1)
%
% with constraint % CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
%
% IMPORTANT!!! FINDORBIT6 is a realistic simulation
% 1. The Frf frequency in the RF cavities (may be different from Frf0)
%    imposes the synchronous condition
%    CT2 - CT1 = C*HarmNumber(1/Frf - 1/Frf0)
% 2. The algorithm numerically calculates
%    6-by-6 Jacobian matrix J6. In order for (J-E) matrix
%    to be non-singular it is NECESSARY to use a realistic
%    PassMethod for cavities with non-zero momentum kick
%    (such as RFCavityPass).
% 3. FINDORBIT6 can find orbits with radiation.
%    In order for the solution to exist the cavity must supply
%    adequate energy compensation.
%    In the simplest case of a single cavity, it must have
%    'Voltage' field set so that Voltage > Erad - energy loss per turn
% 4. FINDORBIT6 starts the search from [ 0 0 0 0 0 0 ]', unless
%    the third argument is specified: FINDORBIT6(RING,REFPTS,GUESS)
%    There exist a family of solutions that correspond to different RF buckets
%    They differ in the 6-th coordinate by C*Nb/Frf. Nb = 1 .. HarmNum-1
% 5. The value of the 6-th coordinate found at the cavity gives
%    the equilibrium RF phase. If there is no radiation the phase is 0;
%
% FINDORBIT6(RING) is 6x1 vector - fixed point at the
%		entrance of the 1-st element of the RING (x,px,y,py,dp,ct)
%
% FINDORBIT6(RING,REFPTS) is 6xLength(REFPTS)
%   array of column vectors - fixed points (x,px,y,py,dp,ct)
%   at the entrance of each element indexed by the REFPTS array.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(RING)+1.
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
% FINDORBIT6(RING,REFPTS,GUESS)
% FINDORBIT6(...,'guess',GUESS)     The search for the fixed point
%	starts from initial condition GUESS. Otherwise the search starts from
%   the synchronous phase. GUESS must be a 6x1 vector.
%
% FINDORBIT6(...,'orbit',ORBIT)     Specify the orbit at the entrance
%   of the ring, if known. FINDORBIT6 will then transfer it to the
%   reference points. ORBIT must be a 6x1 vector.
%
% [ORBIT, FIXEDPOINT] = FINDORBIT6( ... )
%	The optional second return parameter is a 6x1 vector:
%   closed orbit at the entrance of the RING.
%
% See also FINDORBIT4, FINDSYNCORBIT.

if ~iscell(ring)
    error('First argument must be a cell array');
end
[orb6,orbitin] = frequency_control(@xfindorbit6,ring,varargin{:});

    function[orb6,orbitin] = xfindorbit6(ring,varargin)
        [orbitin,varargs]=getoption(varargin,'orbit',[]);
        [refpts,varargs]=getargs(varargs,[],'check',@(arg) isnumeric(arg) || islogical(arg));
        [~,varargs]=getoption(varargs,'is_6d',[]); %% Consume the is_6d option
        if isempty(orbitin)
            orbitin=xorbit_6(ring,varargs{:});
            args={'KeepLattice'};
        else
            args={};
        end

        if islogical(refpts)
            refpts=find(refpts);
        end
        if isempty(refpts)
            % return only the fixed point at the entrance of RING{1}
            orb6=orbitin;
        else
            orb6 = linepass(ring,orbitin,refpts,args{:});
        end
    end
end
