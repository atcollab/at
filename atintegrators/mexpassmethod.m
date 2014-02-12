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
% mex script. All pass-methods #include elempass.h header file
% which uses #if defined(PLATFORM) directive to select
% between platform-specific branches
%
% See also: file:elempass.h

PLATFORMOPTION=[varargin{:}];

%Additional platform-specific options for mex
if ~ispc()
    pdir=fileparts(mfilename('fullpath'));
    switch computer
        case {'SOL2','SOL64'}
            exportarg='-M';
            map='mexFunctionSOL2.map';
            ldflags=['-G -mt ',exportarg,' ' fullfile(pdir,map)];
        case 'GLNX86'
            exportarg='-Wl,--version-script,';
            map='mexFunctionGLNX86.map';
            ldflags=['-pthread -shared -m32 -Wl,--no-undefined ',exportarg,fullfile(pdir,map)];
        case 'GLNXA64'
            exportarg='-Wl,--version-script,';
            map='mexFunctionGLNX86.map';
            ldflags=['-pthread -shared -Wl,--no-undefined ',exportarg,fullfile(pdir,map)];
        case {'MAC','MACI','MACI64'}
            exportarg='-Wl,-exported_symbols_list,';
            map='mexFunctionMAC.map';
    end
    
    if ~exist('verLessThan') || verLessThan('matlab','7.6') %#ok<EXIST>
        PLATFORMOPTION=[PLATFORMOPTION ' LDFLAGS=''',ldflags,''' '];
    elseif verLessThan('matlab','8.3')
        ldf=regexprep(mex.getCompilerConfigurations('C').Details.LinkerFlags,['(' exportarg '\s?)([^\s,]+)'],['$1',fullfile(pdir,map)]);
        PLATFORMOPTION=[PLATFORMOPTION,' LDFLAGS=''',strrep(ldf,'$','\$'),''' '];
    else
        PLATFORMOPTION=[PLATFORMOPTION,' LINKEXPORT=''',exportarg,fullfile(pdir,map),''' '];
    end
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
    if exist(fullfile(pwd,[PM '.c']),'file')
        MEXSTRING = ['mex ',PLATFORMOPTION,PM,'.c'];
        disp(MEXSTRING);
        evalin('base',MEXSTRING);
    else
        disp([PM,'.c',' - NOT FOUND! SKIP']);
    end
end


