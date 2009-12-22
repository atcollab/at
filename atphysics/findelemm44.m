function M44 = findelemm66(ELEM, MethodName, orbit_in);
%FINDELEMM44 numerically finds the 4x4 transfer matrix of an element
%  FINDELEM66(ELEM, METHODNAME, ORBITIN)
%     ELEM          - the element data structure
%     METHODNAME    - name of the pass-method function
%     ORBITIN       - 6-by-1 phase space coordinates at the entrance
%                     The transvese matrix is momentum-dependent,
%                     the 5-th component of ORBITIN is used as the DP value
%
% See also FINDELEMM66

% See if step size for numerical differentiation
% is set globally. Otherwise use 1e-7
global NUMDIFPARAMS
% Transverse
if isfield(NUMDIFPARAMS,'XYStep')
    dt = NUMDIFPARAMS.XYStep';
else
    dt =  1e-7;
end


% Build a diagonal matrix of initial conditions
D4 = [dt*eye(4);zeros(2,4)];
% Add to the orbit_in

RIN = orbit_in*ones(1,8) + [D4, -D4];
% Propagate through the element
ROUT = feval(MethodName,ELEM,RIN);
% Calculate numerical derivative
M44 = [(ROUT(1:4,1:4)-ROUT(1:4,5:8))./(2*dt)];