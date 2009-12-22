function ELEMSEQ = insertindrift(DRIFT0, ELEM1, POS1, varargin)
% INSERTINDRIFT inserts one or more elements into a drift element
%  and returns a sequence (cell array) of elements  ready to to be used 
%  in AT lattice
%
% ELEMSEQ = INSERTELEM(DRIFT0, ELEM1, POS1, ... ELEMN, POSN)
% 
% EXAMPLE: FODO cell
%
% --- 1. Declare elements
%
% D  = atelem('drift','Length',4.5);
% QF = atelem('quad','Length', 1, 'K',  1.234);
% QD = atelem('quad','Length', 1, 'K', -2.345);
%
% --- 2. Insert quadrupoles in the drift;
%
% FODOCELL = insertindrift(D, QF, 0.5, QD, 2, QF, 3.5);
% 
% See also: SPLITELEM


if ~isatelem(DRIFT0)
    error('The first argument must be a valid Accelerator Toolbox drift element');
end

if ~isatelem(ELEM1)
    error('The second argument must be a valid Accelerator Toolbox element to insert');
end

if ~isnumeric(POS1)
    errorstr = sprintf('Incorrect syntax:\n');
    errorstr=[errorstr,sprintf('Elements to inserted must be followed by position [m]\n')];
    errorstr=[errorstr,sprintf('in the argument list: ELEM1, pos1, ... ELEMN, POSN ')];
    error(errorstr);
end
        
    
ELEMSEQ = {};
if POS1>0
    ELEMSEQ{1} = atelem(DRIFT0,'Length',POS1);
elseif POS1<0
    ('Inconsistent lengths and positions cause elements to overlap');
end
    
    
ELEMSEQ{end+1} = ELEM1;
LCUM = POS1+ELEM1.Length;


k = 1;
while k < nargin-3 % Loop to extra
    if isatelem(varargin{k})
        if ~k<nargin & ~isnumeric(varargin{k+1})
            errorstr = sprintf('Incorrect syntax:\n');
            errorstr=[errorstr,sprintf('Elements to inserted must be followed by position [m]\n')];
            errorstr=[errorstr,sprintf('in the argument list: ELEM1, pos1, ... ELEMN, POSN ')];
            error(errorstr);
        else
            
            if (varargin{k+1}-LCUM)>0
                ELEMSEQ{end+1} = atelem(DRIFT0,'Length', varargin{k+1} - LCUM);
            elseif (varargin{k+1}-LCUM)<0
                error('Inconsistent lengths and positions cause elements to overlap');
            end
               
            ELEMSEQ{end+1} = varargin{k};
            LCUM = varargin{k+1}+varargin{k}.Length;
            k = k+2;

            
        end
    else
        error('Incorrect syntax');
    end   
end

if DRIFT0.Length-LCUM > 0
    ELEMSEQ{end+1} = atelem(DRIFT0,'Length',DRIFT0.Length - LCUM);
elseif DRIFT0.Length-LCUM < 0
    error('Inconsistent lengths and positions cause elements to overlap');
end
