function [Rout, varargout] = ringpass(ring, Rin, varargin)
%RINGPASS tracks particles through each element of the cell array RING
% calling the element-specific tracking function specified in the
% RING{i}.PassMethod field.
%
% ROUT=RINGPASS(RING,RIN,NTURNS) tracks particle(s) with initial
%    condition(s) RIN for NTURNS turns
%
%   RING        AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   NTURNS      Number of turns to perform
%   ROUT        6x(N*NTURNS) matrix: output coordinates of N particles at
%               the exit of NTURNS turns
%
% [ROUT, LOST]=RINGPASS(...)
%  Return additionally an information on lost particles
%    LOST	1xN logical vector, indicating lost particles
%    If only one output is given, loss information is saved in
%    global variable LOSSFLAG
%
% [ROUT, LOST, NTURNS]=RINGPASS(...)
%  Return additionally the number of turns performed by each particle
%	NTURNS	1xN vector, number of turns performed
%
% [ROUT, LOSS, NTURNS, LOSSINFO]=RINGPASS(...)
%  Return additional information on lost particles
%   LOSSINFO	1x1 structure with the following fields:
%               turn        1xN vector, turn number where the particle is lost
%               element     1xN vector, element number where the particle is lost
%               coordinate  6xN matrix, coordinates when the particle was
%                           marked as lost
%
% ROUT=RINGPASS(RING,RIN) defaults NTURNS to 1
%
% ROUT=RINGPASS(...,'reuse') with 'reuse' flag is more efficient because
%    it reuses some of the data  and functions stored in the persistent
%    memory from previous calls to RINGPASS.
%
%    !!! In order to use this option, RINGPASS or LINEPASS must first be
%    called without the reuse flag. This will create persistent data structures
%    and keep pointers to pass-method functions.
%
%    !!! RINGPASS(...'reuse') assumes that the number of
%    elements in LINE and pass methods specified in the
%    PassMethod field of each element DO NOT CHANGE between
%    calls. Otherwise, RINGPASS without 'reuse' must be called again.
%    The values of elements fields such as 'Length' or 'K' are allowed to change
%
% ROUT=RINGPASS(...,PREFUNC)
% ROUT=RINGPASS(...,PREFUNC,POSTFUNC)
% ROUT=RINGPASS(...,function_handle.empty,POSTFUNC)
%   PREFUNC and POSTFUNC are function handles, PREFUNC is called
%   immediately before tracking each element, POSTFUNC is called
%   immediately after each element. Functions are called as:
%
%       ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
%
%   and is allowed to modify the particle coordinates
%
% See also: LINEPASS

% Check input arguments
if size(Rin,1)~=6
    error('Matrix of initial conditions, the second argument, must have 6 rows');
end

reuseargs = strcmpi(varargin,'reuse');
newlattice = double(~any(reuseargs));

numericargs = cellfun(@isnumeric,varargin);
nt=find(numericargs,1);
if isempty(nt)
    nturns = 1;
else
    nturns = varargin{nt};
end

funcargs=cellfun(@(arg) isa(arg,'function_handle')||ischar(arg), varargin) & ~reuseargs;

try
    [Rout,lossinfo] = atpass(ring,Rin,newlattice,nturns,length(ring)+1,...
        varargin{funcargs});
    
    if nargout>1;
        if nargout>3, varargout{3}=lossinfo; end
        if nargout>2, varargout{2} = lossinfo.turn; end
        varargout{1} = isfinite(lossinfo.turn);
    else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
        evalin('base','clear LOSSFLAG');
        evalin('base','global LOSSFLAG');
        assignin('base','LOSSFLAG',isfinite(lossinfo.turn));
    end
catch
    error('Atpass:obsolete',['ringpass is now expecting 2 output arguments from atpass.\n',...
        'You may need to call "atmexall" to install the new version']);
end
end
