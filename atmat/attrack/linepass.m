function [Rout,varargout] = linepass(line,Rin,varargin)
%LINEPASS tracks particles through each element of the cell array LINE
% calling the element-specific tracking function specified in the
% LINE{i}.PassMethod field.
%
% ROUT=LINEPASS(LINE,RIN) tracks particle(s) with initial
%    condition(s) RIN for NTURNS turns to the end of the LINE
%
%   LINE        AT lattice
%   RIN         6xN matrix: input coordinates of N particles
%   ROUT        6xN matrix: output coordinates of N particles at
%               the end of LINE
%
% ROUT=LINEPASS(LINE,RIN,REFPTS) also returns intermediate results
%     at the entrance of each element specified in the REFPTS
%
%    REFPTS is an array of increasing indexes that selects elements
%     between 1 and length(LINE)+1.
%     See further explanation of REFPTS in the 'help' for FINDSPOS
%   ROUT        6x(N*length(REFPTS)) matrix: output coordinates of N particles at
%               each reference point
%
%     NOTE:
%     LINEPASS(LINE,RIN,length(LINE)+1) is the same as LINEPASS(LINE,RIN)
%     since the reference point length(LINE)+1 is the exit of the last element
%     LINEPASS(LINE,RIN,1) is a copy of RIN since the
%     reference point 1 is the entrance of the first element
%
% [ROUT, LOST]=LINEPASS(...)
%  Return additionally an information on lost particles
%    LOST	1xN logical vector, indicating lost particles
%    If only one output is given, loss information is saved in
%    global variable LOSSFLAG
%
% [ROUT, LOSS, LOSSINFO]=LINEPASS(...)
%  Return additional information on lost particles
%   LOSSINFO	1x1 structure with the following fields:
%               turn        1xN vector, turn number where the particle is lost
%               element     1xN vector, element number where the particle is lost
%               coordinate  6xN matrix, coordinates when the particle was
%                           marked as lost
%
% ROUT=LINEPASS(...,'reuse') with 'reuse' flag is more efficient because
%    it reuses some of the data  and functions stored in the persistent
%    memory from previous calls to RINGPASS.
%
%    !!! In order to use this option, RINGPASS or LINEPASS must first be
%    called without the reuse flag. This will create persistent data structures
%    and keep pointers to pass-method functions.
%
%    !!! LINEPASS(...'reuse') assumes that the number of
%    elements in LINE and pass methods specified in the
%    PassMethod field of each element DO NOT CHANGE between
%    calls. Otherwise, LINEPASS without 'reuse' must be called again.
%    The values of elements fields such as 'Length' or 'K' are allowed to change
%
% Rfin=LINEPASS(...,PREFUNC)
% Rfin=LINEPASS(...,PREFUNC,POSTFUNC)
% Rfin=LINEPASS(...,function_handle.empty,POSTFUNC)
%    PREFUNC and POSTFUNC are function handles, PREFUNC is called
%    immediately before tracking each element, POSTFUNC is called
%    immediately after each element. Functions are called as:
%
%       ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
%
%   and is allowed to modify the particle coordinates
%
% See also: RINGPASS

% Check input arguments
if size(Rin,1)~=6
    error('Matrix of initial conditions, the second argument, must have 6 rows');
end

reuseargs = strcmpi(varargin,'reuse');
newlattice = double(~any(reuseargs));

numericargs = cellfun(@(arg) isnumeric(arg)||islogical(arg),varargin);
nt=find(numericargs,1);
if isempty(nt)
    refpts = length(line)+1;
elseif islogical(varargin{nt})
    refpts = find(varargin{nt});
else
    refpts = varargin{nt};
end

funcargs=cellfun(@(arg) isa(arg,'function_handle')||ischar(arg), varargin) & ~reuseargs;

try
    [Rout,lossinfo] = atpass(line,Rin,newlattice,1,refpts,varargin{funcargs});
    
    if nargout>1;
        if nargout>2, varargout{2}=lossinfo; end
        varargout{1} = lossinfo.lost;
    else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
        evalin('base','clear LOSSFLAG');
        evalin('base','global LOSSFLAG');
        assignin('base','LOSSFLAG',lossinfo.lost);
    end
catch
    error('Atpass:obsolete',['linepass is now expecting 2 output arguments from atpass.\n',...
        'You may need to call "atmexall" to install the new version']);
end
end
