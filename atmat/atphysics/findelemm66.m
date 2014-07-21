function M66 = findelemm66(ELEM, MethodName, R0)
%FINDELEMM66 numerically finds the 6x6 transfer matrix of an element
%  FINDELEM66(ELEM, METHODNAME, ORBITIN)
%     ELEM          - the element data structure
%     METHODNAME    - name of the pass-method function
%                   (default:  ELEM.PassMethod)
%     ORBITIN       - 6-by-1 phase space coordinates at the entrance
%                   (default: zeros(6,1))
% 
% See also FINDELEMM44

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
% Longitudinal
if isfield(NUMDIFPARAMS,'DPStep')
    dl = NUMDIFPARAMS.DPStep';
else
    % optimal differentiation step - Numerical Recipes
    dl =  6.055454452393343e-006;
end

% Build a diagonal matrix of initial conditions
D6 = [0.5*dt*eye(4),zeros(4,2);zeros(2,4),0.5*dl*eye(2)];
% Add to the orbit_in
RIN = R0(:,ones(1,12)) + [D6, -D6];
% Propagate through the element
ROUT = feval(MethodName,ELEM,RIN);
% Calculate numerical derivative
M66 = [(ROUT(:,1:4)-ROUT(:,7:10))./dt, (ROUT(:,5:6)-ROUT(:,11:12))./dl];
