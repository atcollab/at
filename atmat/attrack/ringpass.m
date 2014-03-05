function [Rout, varargout] = ringpass(ring, Rin, varargin)
%RINGPASS numerically tracks particles through each element in the
% cell array RING calling the element-specific tracking function
% specified in the RING{i}.PassMethod field.
%
% ROUT=RINGPASS(RING,RIN,NUMTURNS) tracks particle(s) with initial
%    condition(s) RIN for NUMTURNS turns
%    Rin is a 6-component  column vector or a 6-by-N matrix
%    made of N columns of different initial conditions.
%    Rfin is 6-by-(number of columns in Rin)*NUMTURNS matrix.
%
% [ROUT, LOSS] =RINGPASS(RING,RIN,NUMTURNS)
%    LOSS is 1 if particle is lost
%    If only one output is given, loss information is saved in
%    global variable LOSSFLAG
%
% ROUT=RINGPASS(RING,RIN) defaults NUMTURNS to 1
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
%    PREFUNC and POSTFUNC are function handles, PREFUNC is called
%    immediately before tracking each element, POSTFUNC is called
%    immediately after each element. Functions are called as:
%       ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
%
% See also: LINEPASS

% Check input arguments
if size(Rin,1)~=6
    error('Matrix of initial conditions, the second argument, must have 6 rows');
end

reuseargs = strcmpi(varargin,'reuse');
if any(reuseargs)
    newlattice = 0;
else
    newlattice = 1;
end

numericargs = cellfun(@isnumeric,varargin);
nt=find(numericargs,1);
if isempty(nt)
    nturns = 1;
else
    nturns = varargin{nt};
end

funcargs=cellfun(@(arg) isa(arg,'function_handle')||ischar(arg), varargin) & ~reuseargs;

[Rout,lossinfo] = atpass(ring,Rin,newlattice,nturns,length(ring)+1,...
    varargin{funcargs});

if nargout>1;
    varargout{1} = isfinite(lossinfo.turn);
    if nargout>2
        varargout{2} = lossinfo.turn;
    end
    if nargout>3
        varargout{3}=lossinfo;
    end
else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
    evalin('base','clear LOSSFLAG');
    evalin('base','global LOSSFLAG');
    assignin('base','LOSSFLAG',isfinite(lossinfo.turn));
end
end
