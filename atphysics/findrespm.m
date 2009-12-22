function C = findrespm(RING, OBSINDEX, PERTURB, PVALUE, varargin)
%FINDRESPM computes the change in the closed orbit due to parameter perturbations
% Two calling syntax options 
% 1. FINDRESPM(RING, OBSINDEX, PERTURBINDEX, PVALUE, 'FIELD', M, N, ORBITFUNCTION, ARGS)
% 2. !!! not implemented yet FINDRESPM(RING, OBSINDEX, PERTURBGROUP, PVALUE, ORBITFUNCTION, ARGS)
%
% RING      - ring lattice
% OBSINDEX  - indexes of elements where the orbit is observed (at the entrance)
% PERTURBINDEX  - Integer indexes of elements whose parameters are perturbed
%                 used with syntax 1 only. 
%             
% PERTURBGROUP  - cell array of AT paramgroups. See ATPARAMGROUP
%               used with syntax 2 only
%
% PVALUE    - amount of peturbation 
%             (Numeric array or scalar if all perturbations are the same magnitude) 
% 
% FIELD,M,N are only use with syntax 1. 
%
% FIELD     - field name of the parameter to perturb (string)
%
% M,N       - index in the matrix, if the field is a matrix
%             For example to perturb the quadrupole field in a
%             multipole element
%             FIELD = 'PolynomB', M = 1, N = 2
%
% ORBITFUNCTION  - specifies which of the FINDORBIT functions is used
%             
%             'findorbit4' (default)
%             'findsyncorbit'
%             'findorbit6'
%             
% ARGS - additioanl arguments may be passsed to some of the FINDORBIT functions
%             findorbit4     - constant momentum error dP
%             findsyncorbit  - fixed orbit lengthening dCT
%              
%
% Returns a 1-by-4 cell array of O-by-P matrixes 
% where O = length(OBSINDEX) and P = length(PERTURB)
% one for each of the close orbit components: X, PX, Y, PY
% See also ATPARAMGROUP, FINDORBIT, FINDORBIT4, FINDORBIT6, FINDSYNCORBIT


O = length(OBSINDEX);
P = length(PERTURB);
C = {zeros(O,P),zeros(O,P),zeros(O,P),zeros(O,P)};

if length(PVALUE) ~= P
    PVALUE = PVALUE(ones(1,P(1)));
end
   

if isnumeric(PERTURB)   % syntax option 1
                        % Integer indexes of perturbed elements. 
                        % More fields must be supplied.
                        % setfield will be used to make perturbations 
if nargin < 7
    error('Incorrect number of inputs');
end

if ~ischar(varargin{1}) % Check that the FIELD argument is a string
    error('The 5-th argument FIELD must be a string');
end
    
if ~isnumeric(varargin{2}) | length(varargin{2})>1 % Check that the M argument is a scalar
    error('The 6-th argument FIELD must be a scalar');
end
M = varargin{2}(1);

if ~isnumeric(varargin{3}) | length(varargin{3})>1 % Check that the M argument is a scalar
    error('The 7-th argument FIELD must be a scalar');
end
N = varargin{3}(1);

if nargin > 7
    ORBITFUNCTION = varargin{4};
else
    ORBITFUNCTION = 'findorbit4';
end
    
   
switch ORBITFUNCTION
case 'findorbit4' 
    orbit_function_handle = @findorbit4;
    if nargin == 9
        orbit_function_args   = {varargin{5}, OBSINDEX};
    else
        orbit_function_args   = {0, OBSINDEX};
    end
case 'findsyncorbit'
    orbit_function_handle = @findsyncorbit;
    if nargin == 9
        orbit_function_args   = {varargin{5}, OBSINDEX};
    else
        orbit_function_args   = {0, OBSINDEX};
    end
case 'findorbit6'
    orbit_function_handle = @findorbit6;
    orbit_function_args   = {OBSINDEX};
otherwise 
    error(['Unknown FINDORBIT function: ',ORBITFUNCTION]);
end

%ORBIT = findorbit4(RING,0,OBSINDEX);


ORBIT = feval(orbit_function_handle,RING,orbit_function_args{:});

mn = {M,N};
for i = 1:P
    oldvalue = getfield(RING{PERTURB(i)},varargin{1},mn);
    RING{PERTURB(i)} = setfield(RING{PERTURB(i)},varargin{1},mn,oldvalue+PVALUE(i));
    ORBITPLUS  = feval(orbit_function_handle,RING,orbit_function_args{:});
    RING{PERTURB(i)} = setfield(RING{PERTURB(i)},varargin{1},mn,oldvalue);
    DORBIT = (ORBITPLUS - ORBIT);
    C{1}(:,i) = DORBIT(1,:);
    C{2}(:,i) = DORBIT(2,:);
    C{3}(:,i) = DORBIT(3,:);
    C{4}(:,i) = DORBIT(4,:);
end
end