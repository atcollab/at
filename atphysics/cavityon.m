function varargout = cavityon(varargin)
%CAVITYON turns Cavities ON
%
% CAVITYON looks for elements that have field Frequency
%    and sets PassMethod for them to ThinCavityPass
% CAVITYON(ENERGY)
%    In addition sets the E0 field of the global variable GLOBVAL
%    to energy - design energy [eV]
%    If GLOBVAL does not exist CAVITYON creates it
%
% See also CAVITYOFF, RADIATIONON, RADIATIONOFF

if ~evalin('base','exist(''THERING'')') | ~evalin('base','isglobal(THERING)');
   error('Global variable THERING could not be found');
end

global THERING

global GLOBVAL
if ~evalin('base','exist(''GLOBVAL'')') | ~evalin('base','isglobal(GLOBVAL)');
    if nargin==0
        error('global variable GLOBVAL does not exist. To create it, use CAVITYON(ENERGY)')
    else
        evalin('base','global GLOBVAL');
        GLOBVAL.E0 = varargin{1};
        disp('Global variable GLOBVAL was created');
    end
else
    if nargin==1
        GLOBVAL.E0 = varargin{1};
    end
    
end

localcavindex = findcells(THERING,'Frequency');

if isempty(localcavindex)
   error('No cavities were found in the lattice');
end

THERING = setcellstruct(THERING,'PassMethod',localcavindex, 'CavityPass');

disp(strcat('Cavities located at index  [',num2str(localcavindex),  ']  were turned ON'))     
if nargout
    varargout{1}=localcavindex;
end