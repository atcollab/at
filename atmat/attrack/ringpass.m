function [Rout, varargout] = ringpass(ring, Rin, nturns, varargin)
%RINGPASS tracks particles through each element of the cell array RING
% calling the element-specific tracking function specified in the
% RING{i}.PassMethod field.
%
% ROUT=RINGPASS(RING,RIN,NTURNS) tracks particle(s) with initial
%    condition(s) RIN for NTURNS turns
%
%   RING        AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   NTURNS      Number of turns to perform (default: 1)
%
%   ROUT        6x(N*NTURNS) matrix: output coordinates of N particles at
%               the exit of each turn
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
% [ROUT, LOSS, NTURNS, LOSSINFO]=RINGPASS(...,'nhist',NHIST,...)
%  Return additional information on lost particles
%   NHIST       number elements before the loss to be traced (default: 1)
%   LOSSINFO	1x1 structure with the following fields:
%               lost                 1xN logical vector, indicating lost particles
%               turn                 1xN vector, turn number where the particle is lost
%               element              1xN vector, element number where the particle is lost
%               coordinates_at_loss  6xN array, coordinates at the exit of
%                                    the element where the particle is lost
%                                    (sixth coordinate is inf if particle is lost in a physical aperture)
%               coordinates          6xNxNHIST array, coordinates at the entrance of the
%                                    LHIST elements before the loss
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
% ROUT=RINGPASS(...,'silent') does not output the particle coordinates at
%    each turn but only at the end of the tracking
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
%   and are allowed to modify the particle coordinates
%
% See also: LINEPASS

% Check input arguments
if size(Rin,1)~=6
    error('Matrix of initial conditions, the second argument, must have 6 rows');
end

if nargin < 3
    nturns = 1;
end
[reuse,args]=getflag(varargin, 'reuse');
[silent,args]=getflag(args, 'silent');
funcargs=cellfun(@(arg) isa(arg,'function_handle'), args);
nhist=getoption(struct(args{~funcargs}), 'nhist',1);

newlattice = double(~reuse);

if silent
    refpts=[];
else
    refpts=length(ring)+1;
end

[prefunc,postfunc]=parseargs({function_handle.empty,function_handle.empty},...
    args(funcargs));

try
    [Rout,lossinfo] = atpass(ring,Rin,newlattice,nturns,refpts,prefunc,postfunc,nhist);
    
    if nargout>1;
        if nargout>3, varargout{3}=lossinfo; end
        if nargout>2, varargout{2} = lossinfo.turn; end
        varargout{1} = lossinfo.lost;
    else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
        evalin('base','clear LOSSFLAG');
        evalin('base','global LOSSFLAG');
        assignin('base','LOSSFLAG',lossinfo.lost);
    end
catch err
    if strcmp(err.identifier,'MATLAB:unassignedOutputs')
        error('Atpass:obsolete',['ringpass is now expecting 2 output arguments from atpass.\n',...
        'You may need to call "atmexall" to install the new version']);
    else
        rethrow(err)
    end
end
end
