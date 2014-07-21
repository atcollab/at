function [M66, varargout] = findm66(LATTICE, varargin)
%FINDM66 numerically finds the 6x6 transfer matrix of an accelerator lattice
%  by differentiation of LINEPASS near the closed orbit
%  FINDM66 uses FINDORBIT6 to search for the closed orbit in 6-d
%  In order for this to work the ring MUST have a CAVITY element
%
% M66 = FINDM66(RING)finds full one-turn 6-by-6
%    matrix at the entrance of the first element
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
%    When REFPTS= [ 1 2 .. ] the fist point is the entrance of the first element
%    and T(:,:,1) - identity matrix
%    When REFPTS= [  .. length(RING)+1] the last point is the exit of the last element
%    and the entrance of the first element after 1 turn: T(:,:, ) = M
%
% [M66,T] = FINDM66(RING,REFPTS,ORBITIN) - Does not search for closed orbit.
%   This syntax is useful to avoid recomputing the closed orbit if it is
%   already known.
%
% [M66,T,orbit] = FINDM66(RING, REFPTS) in addition returns the closed orbit
%    found in the process of lenearization
%
% See also FINDM44, FINDORBIT6

if ~iscell(LATTICE)
    error('First argument must be a cell array');
end

NE = length(LATTICE);

if nargin >= 3 && ~isempty(varargin{2}) % FINDM66(RING,REFPTS,ORBITIN)
    R0 = varargin{2};
else
    R0 = findorbit6(LATTICE);
end
if nargin >= 2  % FINDM66(RING,REFPTS)
    if islogical(varargin{1})
        REFPTS=varargin{1};
        REFPTS(end+1:NE+1)=false;
    elseif isnumeric(varargin{1})
        REFPTS=setelems(false(1,NE+1),varargin{1});
    else
        error('REFPTS must be numeric or logical');
    end
else
    REFPTS=false(1,NE+1);
end
refs=setelems(REFPTS,NE+1);
reqs=REFPTS(refs);

% Determine step size to use for numerical differentiation
global NUMDIFPARAMS
% Transverse
if isfield(NUMDIFPARAMS,'XYStep')
    dt = NUMDIFPARAMS.XYStep';
else
    % optimal differentiation step - Numerical Recipes
    dt =  6.055454452393343e-006;
end
% Longitudinal
if isfield(NUMDIFPARAMS,'DPStep')
    dl = NUMDIFPARAMS.DPStep';
else
    % optimal differentiation step - Numerical Recipes
    dl =  6.055454452393343e-006;
end

% Build a diagonal matrix of initial conditions
D6 = [0.5*dt*eye(4),zeros(4,2);zeros(2,4),0.5*dl*eye(2)];
% Add to the orbit_in. First 12 columns for derivative
% 13-th column is for closed orbit
RIN = R0(:,ones(1,13)) + [D6 -D6 zeros(6,1)];
ROUT = linepass(LATTICE,RIN,refs);
TMAT3 = reshape(ROUT,6,13,[]);
M66 = [(TMAT3(:,1:4,end)-TMAT3(:,7:10,end))./dt,...
    (TMAT3(:,5:6,end)-TMAT3(:,11:12,end))./dl];

if nargout >= 2 % Calculate matrixes at all REFPTS.
    varargout{1} = [(TMAT3(:,1:4,reqs)-TMAT3(:,7:10,reqs))./dt,...
        (TMAT3(:,5:6,reqs)-TMAT3(:,11:12,reqs))./dl];
    % Return closed orbit if requested
    if nargout >= 3
        varargout{2}=squeeze(TMAT3(:,13,reqs));
    end
end

    function mask=setelems(mask,idx)
        mask(idx)=true;
    end

end
