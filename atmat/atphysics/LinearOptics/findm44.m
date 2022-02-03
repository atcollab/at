function [M44, varargout]  = findm44(LATTICE,dp,varargin)
%FINDM44 numerically finds the 4x4 transfer matrix of an accelerator lattice
% for a particle with relative momentum deviation DP
%
% IMPORTANT!!! FINDM44 assumes constant momentum deviation.
%   PassMethod used for any element in the LATTICE SHOULD NOT
%   1.change the longitudinal momentum dP
%     (cavities , magnets with radiation, ...)
%   2.have any time dependence (localized impedance, fast kickers, ...)
%
% M44 = FINDM44(LATTICE,DP) finds the full one-turn
%    matrix at the entrance of the first element
%    !!! With this syntax FINDM44 assumes that the LATTICE
%    is a ring and first finds the closed orbit
%
% [M44,T] = FINDM44(LATTICE,DP,REFPTS) also returns
%    4x4 transfer matrixes  between entrance of
%    the first element and each element indexed by REFPTS.
%    T is 4x4xlength(REFPTS) 3 dimensional array
%    so that the set of indexes (:,:,i) selects the 4-by-4
%    matrix at the i-th reference point.
%
%    Note: REFPTS is an array of increasing indexes that
%    select elements from range 1 to length(LATTICE)+1.
%    See further explanation of REFPTS in the 'help' for FINDSPOS
%    When REFPTS= [ 1 2 .. ] the fist point is the entrance of the
%    first element and T(:,:,1) - identity matrix
%
%    Note: REFPTS is allowed to go 1 point beyond the
%    number of elements. In this case the last point is
%    the EXIT of the last element. If LATTICE is a RING
%    it is also the entrance of the first element
%    after 1 turn: T(:,:,end) = M
%
% [...] = FINDM44(...,'orbit',ORBITIN)
% [...] = FINDM44(RING,DP,REFPTS,ORBITIN)  (Deprecated syntax)
%   Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
%   of initial conditions is used: [x0; px0; y0; py0; DP; 0].
%   The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit if RING is not a
%   ring or to avoid recomputing the closed orbit if is already known.
%
% [...] = FINDM44(...,'full')
%   Same as above except that matrices returned in T are full 1-turn
%   matrices at the entrance of each element indexed by REFPTS.
%
% [...] = FINDM44(...,'ct',CT)
%   Instead of computing the linear optics for  the specified DP/P,
%   computes for the path lenghening specified by CT.
%   The DP argument is ignored.
%
% [M44,T,orbit] = FINDM44(...)
%   In addition returns the closed orbit at the entrance of each element
%   indexed by REFPTS.
%
% See also FINDM66, FINDORBIT4

% *************************************************************************
%   The numerical differentiation in FINDM44 uses symmetric form
%
%         F(x+delta) - F(x-delta)
%       --------------------------------------
%              2*delta
%
%    with optimal differentiation step delta given by !!!! DO LATER
%    The relative error in the derivative computed this way
%    is !!!!!!!!!!!!!!!!! DO LATER
%    Reference: Numerical Recipes.

if ~iscell(LATTICE)
    error('First argument must be a cell array');
end
NE = length(LATTICE);
[fullflag,varargs]=getflag(varargin,'full');
[XYStep,varargs]=getoption(varargs,'XYStep');
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[ct,varargs]=getoption(varargs,'ct',NaN);
[refpts,orbitin,varargs]=getargs(varargs,[],orbitin,'check',@(x) ~(ischar(x) || isstring(x))); %#ok<ASGLU> 

if islogical(refpts)
    refpts(end+1:NE+1)=false;
elseif isnumeric(refpts)
    refpts=setelems(false(1,NE+1),refpts);
else
    error('REFPTS must be numeric or logical');
end

if ~isempty(orbitin)
    if length(orbitin) >= 5
        dp=orbitin(5);
    end
    orbitin = [orbitin(1:4);dp;0];
elseif isnan(ct)
    [~,orbitin] = findorbit4(LATTICE,dp);
else
    [~,orbitin]=findsyncorbit(LATTICE,ct);
end

refs=setelems(refpts,NE+1); % Add end-of-lattice
reqs=refpts(refs);

% Build a diagonal matrix of initial conditions
% scaling=2*XYStep*[1 0.1 1 0.1];
scaling=XYStep*[1 1 1 1];
D4 = [0.5*diag(scaling);zeros(2,4)];
% Add to the orbit_in. First 8 columns for derivative
% 9-th column is for closed orbit
RIN = orbitin + [D4 -D4 zeros(6,1)];
ROUT = linepass(LATTICE,RIN,refs);
TMAT3 = reshape(ROUT(1:4,:),4,9,[]);
M44 = (TMAT3(:,1:4,end)-TMAT3(:,5:8,end))./scaling;

if nargout >= 2 % Calculate matrices at all REFPTS.
    MSTACK = (TMAT3(:,1:4,reqs)-TMAT3(:,5:8,reqs))./scaling;
    
    if fullflag
        S2 = [0 1;-1 0];
        S4 = [S2, zeros(2);zeros(2),S2]; % symplectic identity matrix
        v=cellfun(@rotate,num2cell(MSTACK,[1 2]),'UniformOutput',false);
        varargout{1}=cat(3,v{:});
    else
        varargout{1}=MSTACK;
    end
    % return the closed orbit if requested
    if nargout == 3
        varargout{2}=squeeze(TMAT3(:,9,reqs));
    end
    
end

    function mask=setelems(mask,idx)
        mask(idx)=true;
    end

    function w=rotate(v)
        w=v*M44*S4'*v'*S4;
    end
end
