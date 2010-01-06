function elemstruct = atelem(ELEM,varargin)
%ATELEM makes a new AT element structure from another element, 
% standard AT type, or a template on file. Existing fields are overridden
% with new values. New fields are added and initialized.
%
% NEWELEM = ATELEM(ELEM,'Field1','Value1','Field2', 'Value2', ...)
%  ELEM can be 1) another element structure
%              2) One of the standard AT types
%                 'drift'
%                 'quadrupole'
%                 'sextupole'
%                  ...
%              3) name of a function that returns a template element structure
% 
% Accelerator Toolbox, Version 1.3, 2004-01-26
if isstruct(ELEM)
    elemstruct = ELEM;
    NARG = nargin - 1; 
    k = 1;
    while k < NARG
        NextField = (varargin{k});
        if ischar(NextField)
            elemstruct.(NextField) = varargin{k+1};
        else
            error('Field and value input arguments must come in pairs');
        end
        k = k+2;
    end
elseif ischar(ELEM)
    switch lower(ELEM(1:4)) % Generic
        case 'drif'
            E.FamName = 'None';
            E.Length = 0;
            E.PassMethod = 'DriftPass';
            
        case 'quad'
            E.FamName = 'None';
            E.Length = 0;
            E.K = 0;
            E.PassMethod = 'QuadLinearPass';
            
        case 'sext'
            E.FamName = 'None';
            E.Length = 0;
            E.PolynomB  = [0 0 0];
            E.PolynomA  = [0 0 0];
            E.MaxOrder = 2;
            E.NumIntSteps = 10;
            E.PassMethod = 'StrMPoleSymplectic4Pass';
        case 'mark'
            E.FamName = 'None';
            E.Length = 0;
            E.PassMethod = 'IdentityPass';
            
        case {'bend','rben','sben'}
            E.FamName = 'None';
            E.Length = 0;
            E.BendingAngle = 0;
            E.EntranceAngle =  0;
            E.ExitAngle = 0;
            E.PassMethod = 'BendLinearPass';
            
        case 'corr'
            E.FamName = 'None';
            E.Length = 0;
            E.KickAngle = [0 0];
            E.PassMethod = 'CorrectorPass';
        otherwise % function that returns an at element structure
            if exist(ELEM)==2
                try
                    E = feval(ELEM);
                catch
                    E = [];
                end
                if ~isstruct(E)
                    error(['Function ','''', ELEM,'''', ' does not return a valid element structure']);
                end
            else
                error(['Function ','''',ELEM,'''', ' is not on MATLAB path']);
            end
    end
    elemstruct = atelem(E,varargin{:});
else
    error('First argument must be an AT element or a string keyword')
end