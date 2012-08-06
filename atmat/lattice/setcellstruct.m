function CELLARRAY = setcellstruct(CELLARRAY,field,index,values,varargin)
%SETCELLSTRUCT sets the field values of MATLAB cell array of structures
%
%   Note that the calling syntax must be in the form of assignment:
%   CELLARRAY = SETCELLSTRUCT(CELLARRAY,...)
%   MATLAB does not modify variables that only appear on the right
%   hand side as arguments.
%
% Numeric data
% ---------------------------------------------------------
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M,N)
%	 Sets (M,N) element equal to VALUE when the field data is
%   a matrix. The assigned VALUE may be
%   1. Scalar numeric value to be written to all CELLARRAY elements
%   2. Numeric array of the same length as INDEX array
%
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE,M) can be used
%   for one dimensional vectors
%
% CELLARRAY = SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is the same as
%   SETCELLSTRUCT(CELLARRAY,'field',index,1,1) if the field data
%   is a scalar
%
% Character array
% --------------------------------------------------------------------
% CELLARRAY SETCELLSTRUCT(CELLARRAY,'field',INDEX,VALUE) is a MATLAB
%   cell array of strings when specified fields contain strings.
%   The assignment VALUE may be
%   1. Character string, 
%   2. Character array with the number of rows matching the number of
%       elements in INDEX array,
%   3. Cell array of strings, with either one element or with the same
%       length as index.
%
% See also GETCELLSTRUCT FINDCELLS

if(~iscell(CELLARRAY) || ~isstruct(CELLARRAY{1}) || isempty(CELLARRAY))
    error('The first argument must be a non-empty cell array of structures')
end
% Chechk if the second argument is a string
if(~ischar(field))
    error('The second argument ''field'' must be a character string')
end

% cell array 'mn' here is used as comma separated list
% to pass as an argument to setfield.
if nargin > 4
    mn=varargin;
else
    mn={1};
end

if ischar(values)
    values=cellstr(values);
end

selcells=CELLARRAY(index);
NV = length(selcells);
if length(values)==1
    values=values(ones(NV,1));
end

if isnumeric(values)
    for I=1:NV
        selcells{I} = setfield(selcells{I},field,mn,values(I));
    end
elseif iscell(values)
    for I=1:NV
        selcells{I}.(field)=values{I};
    end
else
    error('The field data must be numeric or character array');
end

CELLARRAY(index)=selcells;


