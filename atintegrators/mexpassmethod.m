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
    EXPORT=' %s';
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
    elseif verLessThan('matlab','9.1')                      % R2016b
        EXPORT=[' LINKEXPORT=''',exportarg,fullfile(pdir,'%s'),''' '];
    else
        EXPORT=[' LINKEXPORTVER=''',exportarg,fullfile(pdir,'%s'),''' '];
    end
end

if ischar(PASSMETHODS) % one file name - convert to a cell array
    if strcmpi(PASSMETHODS,'all')
        % Find all files matching '*Pass.c' wildcard
        D = dir(fullfile(pdir,'*Pass.c'));
        ok=cellfun(@(nm) nm(1)~='.',{D.name});  % Eliminate invisible files
        PASSMETHODS=cellfun(@(nm) strrep(nm,'.c',''),{D(ok).name},'UniformOutput',false);
        try
            generate_passlist(PASSMETHODS);
        catch err
            fprintf(2,'\nCannot generate the list of passmethods: %s\n\n', err.message);
        end
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
                    MEXSTRING = ['mex ',PLATFORMOPTION,'-outdir ',pdir,sprintf(EXPORT,map1),PM];
                    disp(MEXSTRING);
                    evalin('base',MEXSTRING);
                catch
                    MEXSTRING = ['mex ',PLATFORMOPTION,'-outdir ',pdir,sprintf(EXPORT,map2),PM];
                    disp(MEXSTRING);
                    evalin('base',MEXSTRING);
                end
            else
                MEXSTRING = ['mex ',PLATFORMOPTION,'-outdir ',pdir,sprintf(EXPORT,map1),PM];
                disp(MEXSTRING);
                evalin('base',MEXSTRING);
            end
        catch err
            fprintf(2,'Could not compile %s: %s\n',PM,err.message);
        end
    else
        fprintf(2,'%s not found, skip to next\n',PM);
    end
end

    function generate_passlist(passmethods)
        [fid,msg]=fopen(fullfile(pdir,'passmethodlist.m'),'wt');
        if ~isempty(msg)
            error(msg);
        end
        fprintf(fid,'function passmethodlist\n');
        fprintf(fid,'%%PASSMETHODLIST\tUtility function for MATLAB Compiler\n%%\n');
        fprintf(fid,'%%Since passmethods are loaded at run time, the compiler will not include them\n');
        fprintf(fid,'%%unless this function is included in the list of functions to be compiled.\n\n');
        nbytes=cellfun(@(pass) fprintf(fid,'%s\n',pass),passmethods); %#ok<NASGU>
        fprintf(fid,'\nend\n');
        fclose(fid);
    end
end
