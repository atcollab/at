function M44 = findelemm44(ELEM, varargin)
%FINDELEMM44 numerically finds the 4x4 transfer matrix of an element
%  FINDELEMM44(ELEM, METHODNAME)
%     ELEM          - the element data structure
%     METHODNAME    - name of the pass-method function
%                   (default:  ELEM.PassMethod)
%
%  M66=FINDELEMM44(...,'orbit',ORBITIN)  (Deprecated syntax)
%  M66=FINDELEMM44(ELEM, METHODNAME, ORBITIN)
%     ORBITIN       - 6x1 phase space coordinates at the entrance
%                   (default: zeros(6,1))
%                   The transverse matrix is momentum-dependent,
%                   the 5-th component of ORBITIN is used as the DP value
%
% See also FINDELEMM66

[XYStep,varargs]=getoption(varargin,'XYStep');
[R0,varargs]=getoption(varargs,'orbit',zeros(6,1));
[MethodName,R0]=getargs(varargs,ELEM.PassMethod,R0);

% Build a diagonal matrix of initial conditions
% scaling=2*XYStep*[1 0.1 1 0.1];
scaling=2*XYStep*[1 1 1 1];
D4 = [0.5*diag(scaling);zeros(2,4)];
% Add to the orbit_in
RIN = R0 + [D4 -D4];
% Propagate through the element
ROUT = feval(MethodName,ELEM,RIN);
% Calculate numerical derivative
M44 = ((ROUT(1:4,1:4)-ROUT(1:4,5:8))./scaling);
