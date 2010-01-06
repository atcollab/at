function [Rout, varargout] = ringpass(RING, Rin, varargin);
%RINGPASS numerically tracks particles through each element in the 
% cell array RING calling the element-specific tracking function 
% specified in the RING{i}.PassMethod field.
%			
% Rfin=RINGPASS(RING,Rin,NUMTURNS) tracks particle(s) with initial
%    condition(s) Rin for NUMTURNS turns 
%    Rin is a 6-component  column vector or a 6-by-N matrix
%    made of N columns of different initial conditions.
%    Rfin is 6-by-(number of columns in Rin)*NUMTURNS matrix.
%
% [Rfin, LOSS] =RINGPASS(RING,Rin,NUMTURNS)
%    L0SS is 1 if particle is lost
%    If only one output is given, loss information is saved in 
%    global variable LOSSFLAG
%
% Rfin=RINGPASS(RING,Rin) defaults NUMTURNS to 1
%
% Rfin=RINGPASS(RING,Rin,NUMTURNS,'reuse') with 'reuse' flag
%    is more efficient because it reuses some of the data
%    and functions stored in the persistent memory from
%    previous calls to RINGPASS. 
%
%    !!! In order to use this option, RINGPASS or LINEPASS must first be
%    called without the reuse flag. This will 
%    create persistent data structures and keep pointers 
%    to pass-method functions. 
%
%    !!! RINGPASS(...'reuse') assumes that the number of 
%    elements in LINE and pass methods specified in the
%    PassMethod field of each element DO NOT CHANGE between
%    calls. Otherwise, LINEPASS without 'reuse' must 
%    be called again. The values of elements fields such as 'Length' or
%    'K' are allowed to change
%     
%    
% See also: LINEPASS

% Check input arguments
if size(Rin,1)~=6
    error('Matrix of initial conditions, the second argument, must have 6 rows');
end


test = strcmpi(varargin,'reuse');
if any(test)
    NEWLATTICEFLAG = 0;
else 
    NEWLATTICEFLAG = 1;
end

numericargs = varargin(find(~test));


if length(numericargs) > 0
    NUMTURNS = numericargs{1};
else
    NUMTURNS = 1;
end 


Rout = atpass(RING,Rin,NEWLATTICEFLAG,NUMTURNS);

% Find NaN's - lost particles
r = isnan(reshape(Rout(1,:),size(Rin,2),[]))';

if nargout>1;
    
    varargout{1} = any(r,1);
    if nargout>2 % [R, LOSS] = rinpass(..) - two output syntax
        % Count NaN's - conpute turn # when looss occured
        varargout{2} = NUMTURNS - sum(double(r));
    end
else % if no output arguments - create LOSSFLAG, for backward compatibility with AT 1.2
    evalin('base','clear LOSSFLAG');
    evalin('base','global LOSSFLAG');
    assignin('base','LOSSFLAG',any(r));
end

    