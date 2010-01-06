function elemseq = insertelem(ELEM, varargin)
% INSERTELEM inserts one or more elements inside another element
%  and returns the resulting sequence of elements  as a cell array
%  ready to to be used in AT lattice
%
% ELEMSEQ = INSERTELEM(ELEM0, ELEM1, pos1, ... ELEMN, POSN, 'OPT1', .. 'OPTN')

% 
% Use [] in place of an element to be inserted to only split ELEM0

 
% % Optons control fields of elements that are created as a result of slicing
%  options other than inserting into drift are NOT IMPLEMENTED YET - IN AT 1.3 
% 'inherit',  'FIELDNAMES', values of these fields are duplicated in new elements
% 'slice',    'FILEDNAMES'  values are scaled by L/Lo
% 'remove',    Fields are removed from all INNER elements
% 
% 'removefirst'
% 'removelast'
% 
% 
% % Typical option sets
% 'drift' - is the default option useful to insert elemenst inside a drift space
% 'bend'
% 'misalign'


if ~isatelem(ELEM)
    error('The first argument must be a valid Accelerator Toolbox element');
end

% Parse arguments
ELEMS2INSERT = {};

POSITION = {};
LENGTH = {};

OPTIONSET = struct;

k = 1;
while k < nargin
    if isatelem(varargin{k})
        if ~k<nargin & ~isnumeric(varargin{k+1})
            errorstr = sprintf('Incorrect syntax:\n');
            errorstr=[errorstr,sprintf('Elements to be inserted must be followed by position\n')];
            errorstr=[errorstr,sprintf('in the argument list: ELEM1, pos1, ... ELEMN, POSN ')];
            error(errorstr);
        else
            ELEMS2INSERT{end+1} = varargin{k};
            k = k+2;
        end
    elseif ischar(varargin{k})
        switch lower(varargin{k})
            case 'drift'
                OPTIONSET.inherit = {'FamName','PassMethod'};
                OPTIONSET.slice = {'Length'};
                OPTIONSET.removeinner = {};
                OPTIONSET.removefirst = {};
                OPTIONSET.removelast  = {};
            otherwise
                error('Options other than drift are not yet implemented');
        end
    else
        error('Incorrect syntax');
    end
    
end

if isempty(OPTIONSET)
    % 
    OPTIONSET.inherit = {'FamName','PassMethod'};
    OPTIONSET.slice = {'Length'};
    OPTIONSET.removeinner = {};
    OPTIONSET.removefirst = {};
    OPTIONSET.removelast  = {};
% Check if lenghts and positions are consistent
nelemins = length(ELEMS2INSERT);
nelemnew = 1+2*ELEMS2INSERT;

L0 = ELEM.Length;

sp(1) = 0

for k = 1:nelemins
sposition( = 
