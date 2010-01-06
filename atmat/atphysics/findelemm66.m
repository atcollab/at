function M66 = findelemm66(ELEM, MethodName, orbit_in);
%FINDELEMM66 numerically finds the 6x6 transfer matrix of an element
%  FINDELEM66(ELEM, METHODNAME, ORBITIN)
%     ELEM          - the element data structure
%     METHODNAME    - name of the pass-method function
%     ORBITIN       - 6-by-1 phase space coordinates at the entrance
% 
% See also FINDELEMM44

% See if step size for numerical differentiation
% is set globally. Otherwise use 1e-7
global NUMDIFPARAMS
% Transverse
if isfield(NUMDIFPARAMS,'XYStep')
    dt = NUMDIFPARAMS.XYStep';
else
    dt =  1e-7;
end
% Longitudinal
if isfield(NUMDIFPARAMS,'DPStep')
    dl = NUMDIFPARAMS.DPStep';
else
    dl =  1e-7;
end

% Build a diagonal matrix of initial conditions
D6 = [dt*eye(4),zeros(4,2);zeros(2,4), dl*eye(2)];
% Add to the orbit_in
RIN = orbit_in*ones(1,12) + [D6, -D6];
% Propagate through the element
ROUT = feval(MethodName,ELEM,RIN);
% Calculate numerical derivative
M66 = [(ROUT(:,1:4)-ROUT(:,7:10))./(2*dt), (ROUT(:,5:6)-ROUT(:,11:12))./(2*dl)];