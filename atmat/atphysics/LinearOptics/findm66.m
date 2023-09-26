function varargout = findm66(LATTICE, varargin)
%FINDM66 numerically finds the 6x6 transfer matrix of an accelerator lattice
%  by differentiation of LINEPASS near the closed orbit
%  FINDM66 uses FINDORBIT6 to search for the closed orbit in 6-d
%  In order for this to work the ring MUST have a CAVITY element
%
% M66 = FINDM66(RING) finds the full one-turn 6-by-6
%    matrix at the entrance of the first element
%
%[...]=FINDM66(RING,...,'dp',DP) Specify the momentum deviation when
%   radiation is OFF (default: 0)
%
%[...]=FINDM66(RING,...,'dct',DCT) Specify the path lengthening when
%   radiation is OFF (default: 0)
%
%[...]=FINDM66(RING,...,'df',DF) Specify the RF frequency deviation
%   radiation is OFF (default: 0)
%
%[...]=FINDM66(RING,...,'orbit',ORBIT) Specify the orbit at the entrance
%   of the ring, if known.
%
% [M66,T] = FINDM66(RING,REFPTS) in addition to M finds
%    6-by-6 transfer matrixes  between entrances of
%    the first element and each element indexed by REFPTS.
%    T is 6-by-6-by-length(REFPTS) 3 dimentional array.
%
%    REFPTS is an array of increasing indexes that  select elements
%    from the range 1 to length(RING)+1.
%    See further explanation of REFPTS in the 'help' for FINDSPOS
%
%    Note:
%    When REFPTS= [ 1 2 .. ] the first point is the entrance of the first element
%    and T(:,:,1) - identity matrix
%    When REFPTS= [  .. length(RING)+1] the last point is the exit of the last element
%    and the entrance of the first element after 1 turn: T(:,:, ) = M
%
% [...] = FINDM66(RING,REFPTS,ORBITIN)    (Deprecated syntax)
% [...] = FINDM66(...,'orbit',ORBITIN)
%   Do not search for closed orbit. This syntax is useful to avoid
%   recomputing the closed orbit if it is already known.
%
% [M66,T,orbit] = FINDM66(...)
%   In addition returns the closed orbit at the entrance of each element
%   indexed by REFPTS.
%
%
% See also FINDM44, FINDORBIT6

if ~iscell(LATTICE)
    error('First argument must be a cell array');
end
[varargout{1:nargout}] = frequency_control(@xfindm66,LATTICE,varargin{:});

    function [M66, varargout] = xfindm66(LATTICE, varargin)
        NE = length(LATTICE);
        [XYStep,varargs]=getoption(varargin,'XYStep');	% Step size for numerical differentiation	%1.e-8
        [DPStep,varargs]=getoption(varargs,'DPStep');	% Step size for numerical differentiation	%1.e-6
        [orbitin,varargs]=getoption(varargs,'orbit',[]);
        [dpargs,varargs]=getoption(varargs,{'dp','dct','df'});
        [is_6d,varargs]=getoption(varargs,'is_6d',[]); % Always set by frequency_control
        [refpts,orbitin,varargs]=getargs(varargs,[],orbitin,'check',@(x) ~(ischar(x) || isstring(x))); %#ok<ASGLU>

        if islogical(refpts)
            refpts(end+1:NE+1)=false;
        elseif isnumeric(refpts)
            refpts=setelems(false(1,NE+1),refpts);
        else
            error('REFPTS must be numeric or logical');
        end

        if isempty(orbitin)
            [~, orbitin] = findorbit(LATTICE,'XYStep',XYStep,'DPStep',DPStep,'is_6d',is_6d,dpargs{:});
        end

        refs=setelems(refpts,NE+1);
        reqs=refpts(refs);

        % Build a diagonal matrix of initial conditions
        %scaling=2*XYStep*[1 0.1 1 0.1 1 1];
        scaling=XYStep*[1 1 1 1 0 0] + DPStep*[0 0 0 0 1 1];
        D6 = 0.5*diag(scaling);
        % Add to the orbit_in. First 12 columns for derivative
        % 13-th column is for closed orbit
        RIN = orbitin + [D6 -D6 zeros(6,1)];
        ROUT = linepass(LATTICE,RIN,refs);
        TMAT3 = reshape(ROUT,6,13,[]);
        M66 = (TMAT3(:,1:6,end)-TMAT3(:,7:12,end))./scaling;

        if nargout >= 2 % Calculate matrices at all REFPTS.
            varargout{1} = (TMAT3(:,1:6,reqs)-TMAT3(:,7:12,reqs))./scaling;
            % Return closed orbit if requested
            if nargout >= 3
                varargout{2}=squeeze(TMAT3(:,13,reqs));
            end
        end

        function mask=setelems(mask,idx)
            mask(idx)=true;
        end

    end
end
