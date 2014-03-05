function [Rout,varargout] = linepass(line,Rin,varargin)
%LINEPASS tracks particles through a sequence of elements 
% in the cell array LINE. For each element it calls 
% the pass-method specified in the 'PassMethod' field. 
%		
% Rout = LINEPASS(LINE,Rin) tracks particle(s) to the end 
%     of the LINE. Rin and Rout are 6-by-1 
%     column vectors or 6-by-N matrixes
%     Each column represents a different initial condition or particle.
%
% Rout = LINEPASS(LINE,Rin,REFPTS) also returns intermediate results 
%     at the entrance of each element specified in the REFPTS
%
%     REFPTS is an array of increasing indexes that  selects elements 
%     between 1 and length(LINE)+1. 
%     See further explanation of REFPTS in the 'help' for FINDSPOS
%     
%     NOTE:
%     LINEPASS(LINE,Rin,length(LINE)+1) is the same as  LINEPASS(LINE,Rin)
%     since the reference point length(LINE)+1 is the exit of the last element
%     LINEPASS(LINE,Rin,1) is a copy of Rin since the 
%     reference point 1 is the entrance of the first element
%     
%     OUTPUT FORMAT:
%     Rout is 6-by-(number of columns in Rin)*length(REFPTS) matrix
%     where blocks 6-by-(number of columns in Rin) corresponds 
%     to different REFPTS
%     FOR EXAMPLE:
%     if Rin is 6-by-2 maid of two 6-by-1 column vectors [Rin1, Rin2]
%     and REFPTS = [N1 N2 N3] so that N1<N2<N3
%     the output is [Rout1(N1) Rout2(N1) Rout1(N2) Rout2(N2) Rout1(N3) Rout2(N3)]
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
%       ROUT=FUNC(ELEMENT, RIN, NTURN, NELEMENT)
%     
% See also: RINGPASS FINDSPOS

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
    refpts = length(line)+1;
else
    refpts = varargin{nt};
end

funcargs=cellfun(@(arg) isa(arg,'function_handle')||ischar(arg), varargin) & ~reuseargs;

[Rout,lossinfo] = atpass(line,Rin,newlattice,1,refpts,varargin{funcargs});

if nargout>1;
    varargout{1} = isfinite(lossinfo.turn);
    if nargout>2
        varargout{2}=lossinfo;
    end
else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
    evalin('base','clear LOSSFLAG');
    evalin('base','global LOSSFLAG');
    assignin('base','LOSSFLAG',isfinite(lossinfo.turn));
end
end
