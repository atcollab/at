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
            elemstruct=atdrift('',varargin{:});
        case 'quad'
            elemstruct=atquadrupole('',varargin{:});
        case 'sext'
            elemstruct=atsextupole('',varargin{:});
        case 'mark'
            elemstruct=atmarker('',varargin{:});
        case {'bend','sben'}
            elemstruct=atsbend('',varargin{:});
        case 'rben'
            elemstruct=atrbend('',varargin{:});
        case 'corr'
            elemstruct=atcorrector('',varargin{:});
        otherwise % function that returns an at element structure
            error('AT:atelem:UnknownType',['Unknown element type: ' elem]);
    end
else
    error('First argument must be an AT element or a string keyword')
end
