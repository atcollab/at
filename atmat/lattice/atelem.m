function elemstruct = atelem(elem,varargin)
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
%                 'bend','rbend','sbend'
%                 'marker'
%                 'corrector'
%                  ...

if isstruct(elem)
    elemstruct = elem;
    for k=1:2:length(varargin)
        NextField = (varargin{k});
        if ischar(NextField)
            elemstruct.(NextField) = varargin{k+1};
        else
            error('Field and value input arguments must come in pairs');
        end
    end
elseif ischar(elem)
    switch lower(elem(1:4)) % Generic
        case 'drif'
            E.Length = 0;
            E.PassMethod = 'DriftPass';
            
        case 'quad'
            E.Length = 0;
            E.K = 0;
            E.PassMethod = 'QuadLinearPass';
            
        case 'sext'
            E.Length = 0;
            E.PolynomB  = [0 0 0];
            E.PolynomA  = [0 0 0];
            E.MaxOrder = 2;
            E.NumIntSteps = 10;
            E.PassMethod = 'StrMPoleSymplectic4Pass';
            
        case 'mark'
            E.Length = 0;
            E.PassMethod = 'IdentityPass';
            
        case {'bend','rben','sben'}
            E.Length = 0;
            E.BendingAngle = 0;
            E.EntranceAngle =  0;
            E.ExitAngle = 0;
            E.PassMethod = 'BendLinearPass';
            
        case 'corr'
            E.Length = 0;
            E.KickAngle = [0 0];
            E.PassMethod = 'CorrectorPass';
            
        otherwise % function that returns an at element structure
            error('AT:atelem:UnknownType',['Unknown element type: ' elem]);
    end
    elemstruct = atelem(E,varargin{:});
else
    error('First argument must be an AT element or a string keyword')
end
