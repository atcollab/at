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
%[...]=FINDORBIT(RING,...,'df',DF) Specify the RF frequency deviation
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

[orbs,orbitin]=frequency_control(@do,ring,varargin{:});

    function [orbs,orb0]=do(ring,varargin)
        [orb0,varargs]=getoption(varargin,'orbit',[]);
        [refpts,varargs]=getargs(varargs,[],'check',@(arg) isnumeric(arg) || islogical(arg));
        [dp,varargs]=getoption(varargs,'dp',NaN);
        [dct,varargs]=getoption(varargs,'dct',NaN);
        [df,varargs]=getoption(varargs,'df',NaN);
        [is_6d,varargs]=getoption(varargs,'is_6d',[]); % Always set by frequency_control
        if isempty(orb0)
            if is_6d   % Radiation ON: 6D orbit
                orb0=xorbit_6(ring,varargs{:});
            elseif isfinite(df)
                [cell_l,cell_frev,cell_h]=atGetRingProperties(ring,'cell_length',...
                    'cell_harmnumber','cell_revolution_frequency');
                dct=-cell_l*df/(cell_frev*cell_h+df);
                orb0=xorbit_ct(ring,dct,varargs{:});
            elseif isfinite(dct)
                orb0=xorbit_ct(ring,dct,varargs{:});
            else
                orb0=xorbit_dp(ring,dp,varargs{:});
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
            orbs=orb0;
        else
            orbs=linepass(ring,orb0,refpts,args{:});
        end
    end
end
