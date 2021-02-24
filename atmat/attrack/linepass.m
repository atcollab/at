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
%
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
% [ROUT, LOSS, LOSSINFO]=LINEPASS(...,'nhist',NHIST,...)
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
% ROUT=LINEPASS(...,'KeepLattice') Tracking with the 'KeepLattice' flag is
%   more efficient because it reuses persistent data structures stored in
%   memory in previous calls to LINEPASS.
%
%	!!! In order to use this option, LINEPASS must first be called
%	without the 'KeepLattice' flag. It then assumes that the elements in LINE
% 	DO NOT CHANGE between calls. Otherwise, LINEPASS must be called again
%   without 'KeepLattice'.
%
% ROUT=LINEPASS(...,'reuse') is kept for compatibilty with previous
% versions. It has no effect.
%
% Rfin=LINEPASS(...,PREFUNC)
% Rfin=LINEPASS(...,PREFUNC,POSTFUNC)
% Rfin=LINEPASS(...,cell(0),POSTFUNC)
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
[keeplattice,args]=getflag(varargin, 'KeepLattice');
[dummy,args]=getflag(args,'reuse');	%#ok<ASGLU> % Kept for compatibility and ignored
[nhist,args]=getoption(args,'nhist',1);
[omp_num_threads,args]=getoption(args,'omp_num_threads');
funcargs=cellfun(@(arg) isa(arg,'function_handle'), args);
refpts=getargs(args(~funcargs),length(line)+1);
[prefunc,postfunc]=getargs(args(funcargs),cell(0),cell(0));

if islogical(refpts)
    refpts=find(refpts);
end
newlattice = double(~keeplattice);

try
    [Rout,lossinfo] = atpass(line,Rin,newlattice,1,refpts,prefunc,postfunc,nhist,omp_num_threads);
    
    if nargout>1
        if nargout>2, varargout{2}=lossinfo; end
        varargout{1} = lossinfo.lost;
    else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
        evalin('base','clear LOSSFLAG');
        evalin('base','global LOSSFLAG');
        assignin('base','LOSSFLAG',lossinfo.lost);
    end
catch err
    if strcmp(err.identifier,'MATLAB:unassignedOutputs')
        error('Atpass:obsolete',['linepass is now expecting 2 output arguments from atpass.\n',...
            'You should call "atmexall" to install the new version']);
    else
        rethrow(err)
    end
end
end
