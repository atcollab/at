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
pdir=fileparts(mfilename('fullpath'));
if ispc()
    EXPORT='%s';
    map1='';
else
    switch computer
        case {'SOL2','SOL64'}
            exportarg='-M';
            map1='mexFunctionSOL2.map';
            ldflags=['-G -mt ',exportarg,' '];
        case 'GLNX86'
            exportarg='-Wl,--version-script,';
            map1='mexFunctionGLNX86.map';
            ldflags=['-pthread -shared -m32 -Wl,--no-undefined ',exportarg];
        case 'GLNXA64'
            exportarg='-Wl,--version-script,';
            map1='mexFunctionGLNX86.map';
            ldflags=['-pthread -shared -Wl,--no-undefined ',exportarg];
        case {'MAC','MACI','MACI64'}
            exportarg='-Wl,-exported_symbols_list,';
            map1='mexFunctionMAC.map';
            map2='mexFunctionMAC2.map';
    end
    
    if ~exist('verLessThan') || verLessThan('matlab','7.6') %#ok<EXIST> R2008a
        EXPORT=[' LDFLAGS=''',ldflags,fullfile(pdir,'%s'),''' '];
    elseif verLessThan('matlab','8.3')                      % R2014a
        ldf=regexprep(mex.getCompilerConfigurations('C').Details.LinkerFlags,['(' exportarg '\s?)([^\s,]+)'],['$1',fullfile(pdir,'%s')]);
        EXPORT=[' LDFLAGS=''',strrep(ldf,'$','\\$'),''' '];
    else
        EXPORT=[' LINKEXPORT=''',exportarg,fullfile(pdir,'%s'),''' '];
    end
end

if ischar(PASSMETHODS) % one file name - convert to a cell array
    if strcmpi(PASSMETHODS,'all')
        % Find all files matching '*Pass.c' wildcard
        D = dir(fullfile(pdir,'*Pass.c'));
        ok=cellfun(@(nm) nm(1)~='.',{D.name});  % Eliminate invisible files
        PASSMETHODS=cellfun(@(nm) strrep(nm,'.c',''),{D(ok).name},'UniformOutput',false);
    else % Mex a single specifie pass-method
        PASSMETHODS={PASSMETHODS};
    end
end

for i = 1:length(PASSMETHODS)
    PM = fullfile(pdir,[PASSMETHODS{i} '.c']);
    evalin('base',['clear ',PASSMETHODS{i}]);
    if exist(PM,'file')
        try
            if exist('map2','var')
                try
                    MEXSTRING = ['mex ',PLATFORMOPTION,sprintf(EXPORT,map1),'-outdir ',pdir,' ',PM];
                    disp(MEXSTRING);
                    evalin('base',MEXSTRING);
                catch
                    MEXSTRING = ['mex ',PLATFORMOPTION,sprintf(EXPORT,map2),'-outdir ',pdir,' ',PM];
                    disp(MEXSTRING);
                    evalin('base',MEXSTRING);
                end
            else
                MEXSTRING = ['mex ',PLATFORMOPTION,sprintf(EXPORT,map1),'-outdir ',pdir,' ',PM];
                disp(MEXSTRING);
                evalin('base',MEXSTRING);
            end
        catch err
            disp(['Could not compile ' PM char(10) err.message]);
        end
    else
        disp([PM,'.c',' - NOT FOUND! SKIP']);
    end
end
end
