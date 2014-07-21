function M44 = findelemm44(ELEM, MethodName, R0)
%FINDELEMM44 numerically finds the 4x4 transfer matrix of an element
%  FINDELEM66(ELEM, METHODNAME, ORBITIN)
%     ELEM          - the element data structure
%     METHODNAME    - name of the pass-method function
%                   (default:  ELEM.PassMethod)
%     ORBITIN       - 6-by-1 phase space coordinates at the entrance
%                   (default: zeros(6,1))
%                   The transverse matrix is momentum-dependent,
%                   the 5-th component of ORBITIN is used as the DP value
%
% See also FINDELEMM66

if (nargin < 3) || isempty(R0),  R0 = zeros(6,1); end
if (nargin < 2) || isempty(MethodName),  MethodName=ELEM.PassMethod; end

% Determine step size to use for numerical differentiation
global NUMDIFPARAMS
% Transverse
if isfield(NUMDIFPARAMS,'XYStep')
    dt = NUMDIFPARAMS.XYStep';
else
    % optimal differentiation step - Numerical Recipes
    dt =  6.055454452393343e-006;
end

% Build a diagonal matrix of initial conditions
D4 = [0.5*dt*eye(4);zeros(2,4)];
% Add to the orbit_in
RIN = R0(:,ones(1,8)) + [D4 -D4];
% Propagate through the element
ROUT = feval(MethodName,ELEM,RIN);
% Calculate numerical derivative
M44 = ((ROUT(1:4,1:4)-ROUT(1:4,5:8))./dt);
