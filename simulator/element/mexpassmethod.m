function mexpassmethod(PASSMETHODS, varargin)
%MEXPASSMETHOD builds pass-method mex-files from C files
%
% PASSMETHODS argument can be:
%  Single pass-method name - the same as the C file name without '.c'
%  Cell array of strings containing pass-method names
%  'all' - option automatically detects all C files matching *Pass.c pattern
%
% The second argument is a list of options passed to the 'mex' script
% 
% Examples: mexpassmethod('DriftPass','-v')
%           mexpassmethod('all','-argcheck')
%           mexpassmethod({'DriftPass','BendLinearPass'})
%
% Note:  MEXPASSMETHOD automatically determines the host 
% platform and costructs -D<PLATFORM> option to feed to the 
% mex script. All pass-methods #incude elempass.h header file
% which uses #if defined(PLATFORM) directive to select
% between platform-specific branches
%
% See also: file:elempass.h

PLATFORMOPTION = ['-D',computer,' '];
CURRENTDIR = pwd;
cd(fileparts(which('DriftPass')));

tmpfile = 0;
%Additional platform-specific options for mex
switch computer
case 'SOL2'
    PLATFORMOPTION = [PLATFORMOPTION,'LDFLAGS=''-shared -W1,-M,',atroot,'/simulator/element/mexFunctionSOL2.map''',' '];
case 'GLNX86'
    PLATFORMOPTION = [PLATFORMOPTION,'LDFLAGS=''-pthread -shared -m32 -Wl,--version-script,',atroot,'/simulator/element/mexFunctionGLNX86.map''',' '];  
end


if ischar(PASSMETHODS) % one file name - convert to a cell array
    if strcmpi(PASSMETHODS,'all')
        % Find all files matching '*Pass.c' wildcard   
        D = dir('*Pass.c');
        PASSMETHODS = cell(size(D));
        for i = 1:length(D)
            PASSMETHODS{i} = strrep(D(i).name,'.c','');
        end
    else % Mex a single specifie pass-method
        PASSMETHODS={PASSMETHODS};
    end
end

for i = 1:length(PASSMETHODS)

        PM = PASSMETHODS{i};
        evalin('base',['clear ',PM]);
        MEXSTRING = ['mex ',PLATFORMOPTION];
        if nargin==2
            MEXSTRING = [MEXSTRING,varargin{1},' '];
        end
        MEXSTRING = [MEXSTRING, PM,'.c '];
        
        %message = sprintf('%s\n',MEXSTRING);
        %disp(message);
        
        if exist([pwd,'\',PM,'.c'],'file') | exist ([pwd,'/',PM,'.c'],'file') 
            disp(MEXSTRING);
            evalin('base',MEXSTRING);
        else 
            disp([PM,'.c',' - NOT FOUND! SKIP']); 
        end
        
end

cd(CURRENTDIR);



