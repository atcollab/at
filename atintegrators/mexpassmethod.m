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

% Optional flag to turn on multi-threading feature (openmp) in some of the
% pass methods to decrease the tracking time for multi-particle
% simulations. The if the number of particle is less than
% OMP_PARTICLE_THRESHOLD (defined in atcommon.h) then only a single thread
% is used. The number of threads used can be manually set using the
% environment variable OMP_NUM_THREADS.
% Originally intriduced by Xiabiao Huang (7/12/2010).
opt_parallel = false;

PLATFORMOPTION=[varargin{:}];

%Additional platform-specific options for mex
pdir=fileparts(mfilename('fullpath'));
if ispc()
    EXPORT=' %s';
    map1='';
else
    % Starting from R2016b, Matlab introduced a new entry point in MEX-files
    % The "*.mapext" files define this new entry point
    if ~exist('verLessThan') || verLessThan('matlab','9.1') %#ok<EXIST>
        mapformat=fullfile(pdir,'%s');
    else
        mapformat=fullfile(pdir,'%sext');
    end
    
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
            map1='trackFunctionMAC.map';
            map2='passFunctionMAC.map';
    end
    
    if opt_parallel
        % Add library flags to enable openmp
        exportarg = [' -fopenmp ', exportarg];
    end
    
    if ~exist('verLessThan') || verLessThan('matlab','7.6') %#ok<EXIST> R2008a
        EXPORT=[' LDFLAGS=''',ldflags,mapformat,''' '];
    elseif verLessThan('matlab','8.3')                      % R2014a
        ldf=regexprep(mex.getCompilerConfigurations('C').Details.LinkerFlags,['(' exportarg '\s?)([^\s,]+)'],['$1',mapformat]);
        EXPORT=[' LDFLAGS=''',strrep(ldf,'$','\\$'),''' '];
    elseif verLessThan('matlab','9.1')                      % R2016b
        EXPORT=[' LINKEXPORT=''',exportarg,mapformat,''' '];
    else
        if ismac()  % Correct a bug in Mac setup which uses both LINKEXPORT and LINKEXPORTVER
            PLATFORMOPTION=[PLATFORMOPTION ...
                'LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch x86_64 -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -Wl,-syslibroot,$ISYSROOT'' ' ...
                'CMDLINE200=''$LD $LDFLAGS $LDBUNDLE $LINKOPTIM $LINKEXPORTVER $OBJS $CLIBS $LINKLIBS -o $EXE'' '...
                ];
        end
        EXPORT=[' LINKEXPORTVER=''',exportarg,mapformat,''' '];
    end
end

if opt_parallel
    % Add compiler flags to enable openmp
    if ispc()
        warning('OpenMP available with MSVC or V Studio compilers.');
        EXPORT = [' COMPFLAGS=''$COMPFLAGS /openmp'' ', EXPORT];
    else
        switch computer
            case {'SOL2','SOL64'}
                EXPORT = [' CFLAGS=''$CFLAGS -xopenmp'' ', EXPORT];
            case {'MAC','MACI','MACI64'}
                warning('OpenMP may not be supportded on MAC. *Untested*');
                EXPORT = [' CFLAGS=''$CFLAGS -fopenmp'' ', EXPORT];
            otherwise
                EXPORT = [' CFLAGS=''$CFLAGS -fopenmp'' ', EXPORT];
        end
    end
end
    

if ischar(PASSMETHODS) % one file name - convert to a cell array
    if strcmpi(PASSMETHODS,'all')
        % Find all files matching '*Pass.c' wildcard
        cfiles = dir(fullfile(pdir,'*Pass.c'));
        passmethods = {cfiles.name};
        % Find all files matching '*Pass.cc' wildcard
        ccfiles = dir(fullfile(pdir, '*Pass.cc'));
        passmethods = [passmethods ccfiles.name];
        % Eliminate invisible files
        ok=cellfun(@(nm) nm(1)~='.',passmethods,'UniformOutput',false);
        passmethods = passmethods(cell2mat(ok));
        try
            generate_passlist(pdir,passmethods);
        catch err
            fprintf(2,'\nCannot generate the list of passmethods: %s\n\n', err.message);
        end
        PASSMETHODS=passmethods;
    else % Mex a single specifie pass-method
        PASSMETHODS={PASSMETHODS};
    end
end

for i = 1:length(PASSMETHODS)
    PM = fullfile(pdir,[PASSMETHODS{i}]);
    if exist(PM,'file') || exist(PM2, 'file')
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

    function generate_passlist(pdir,passmethods)
        % Remove trailing '.c' or '.cc' from each passmethod using regexp
        % Literally, replace dot plus any number of cs at the end of a
        % passmethod with an empty string.
        passmethodnames = cellfun(@(pass) regexprep(pass, '\.c*$', ''),passmethods, 'UniformOutput', false);
        [fid,msg]=fopen(fullfile(pdir,'passmethodlist.m'),'wt');
        if ~isempty(msg)
            error(msg);
        end
        fprintf(fid,'function passmethodlist\n');
        fprintf(fid,'%%PASSMETHODLIST\tUtility function for MATLAB Compiler\n%%\n');
        fprintf(fid,'%%Since passmethods are loaded at run time, the compiler will not include them\n');
        fprintf(fid,'%%unless this function is included in the list of functions to be compiled.\n\n');

        nbytes=cellfun(@(pass) fprintf(fid,'%s\n',pass),passmethodnames); %#ok<NASGU>
        fprintf(fid,'\nend\n');
        fclose(fid);
    end
end
