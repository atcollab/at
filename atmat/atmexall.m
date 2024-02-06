function atmexall(varargin)
%ATMEXALL build all AT platform dependent mex-files from C-sources
%
%ATMEXALL option1 ... optionN
%
% AT Options:
%
%	-missing    Build only the outdated components
%   -fail       Throw an exception if compiling any passmethod fails
%               (By defaults compilation goes on)
%	-openmp     Build the integrators for OpenMP parallelisation
%	-cuda CUDA_PATH Build the GPU tracking support using Cuda
%   -c_only     Do no compile C++ passmethods
%   -DOMP_PARTICLE_THRESHOLD=n
%               Set the parallelisation threshold to n particles
%               (Default 10)
%
% Options forwarded to the mex command:
%
%   -v          Verbose output
%   -g          Compile with debug options
%   -O          Optimize the object code (Default)
%   -n          Display the generated command without executing
%   ...
%

pdir=fullfile(fileparts(atroot),'atintegrators');
[openmp,varargs]=getflag(varargin,'-openmp');
[cuda,varargs]=getoption(varargs,'-cuda','None');
[miss_only,varargs]=getflag(varargs,'-missing');
[c_only,varargs]=getflag(varargs,'-c_only');
[fail,varargs]=getflag(varargs,'-fail');
force=~miss_only;

atoptions={['-D',computer]};

if isunix && ~ismac
    LIBDL={'-ldl'};
else
    LIBDL={};
end

ompflags={};
cflags={};
ldflags={};
if openmp
    if ispc()
        cflags={'/openmp'};
        ompflags={};
    elseif ismac()
        ompinc = select_omp();
        cflags={'-Xpreprocessor', '-fopenmp'};
        ompflags={ompinc,...
             sprintf('-L"%s"',fullfile(matlabroot,'sys','os',computer('arch'))),...
            '-liomp5'};
    else
        cflags={'-fopenmp'};
        ompflags={...
            sprintf('-L"%s"',fullfile(matlabroot,'sys','os',computer('arch'))),...
            '-liomp5'};
    end
end

if ispc()
    ompoptions=pc_flags(ompflags, cflags, ldflags);
    map1=pc_flags(ompflags, cflags, ldflags);
elseif ismac()
    ompoptions=unix_flags(ompflags, cflags, ldflags);
    map1=unix_flags(ompflags, cflags, ldflags, '-Wl,-exported_symbols_list,', 'trackFunctionMAC.map');
    map2=unix_flags(ompflags, cflags, ldflags, '-Wl,-exported_symbols_list,', 'passFunctionMAC.map');
else
    ompoptions=unix_flags(ompflags, cflags, ldflags);
    map1=unix_flags(ompflags, cflags, ldflags, '-Wl,--version-script,', 'mexFunctionGLNX86.map');
end

try
    if ~verLessThan('matlab','9.4')         % >= R2018a
        atoptions = [atoptions,{'-R2018a'}];
    elseif ~verLessThan('matlab','7.11')	% >= R2010b
        atoptions = [atoptions,{'-largeArrayDims'}];
    end
    if ~verLessThan('matlab','8.3')         % >= R2014a
        atoptions = [atoptions,{'-silent'}];
    end
catch
end

lastwarn('');
passinclude = ['-I' pdir];
alloptions=[atoptions varargs];

% atpass
cdir=fullfile(atroot,'attrack','');
compile([alloptions, {passinclude}, LIBDL, ompoptions], fullfile(cdir,'atpass.c'));
compile([atoptions, ompoptions],fullfile(cdir,'coptions.c'))

% gpuextensions
if ~strcmp(cuda,'None')
    gpudir=fullfile(fileparts(atroot),'atgpu','');
    generate_gpuheader(gpudir)
    if ispc()
        % TODO
        error('AT:atmexall', 'GPU windows not supported');
    elseif ismac()
        % TODO
        error('AT:atmexall', 'GPU ismac not supported');
    else
        gpuflags = {sprintf('-I"%s"',gpudir),...
                    sprintf('-I"%s/include"',cuda),...
                    sprintf('-L"%s/lib64"',cuda),...
                    sprintf('LDFLAGS=$LDFLAGS -Wl,-rpath,"%s/lib64"',cuda),...
                    '-DCUDA'};
    end
    compile([alloptions, {passinclude}, gpuflags], ...
        fullfile(cdir,'gpuinfo.cpp'),...
        fullfile(gpudir,'MatlabInterface.cpp'), ...
        fullfile(gpudir,'AbstractInterface.cpp'), ...
        fullfile(gpudir,'CudaGPU.cpp'), ...
        fullfile(gpudir,'AbstractGPU.cpp'), ...
        '-lcuda','-lnvrtc');
    compile([alloptions, {passinclude}, gpuflags], ...
        fullfile(cdir,'gpupass.cpp'),...
        fullfile(gpudir,'AbstractGPU.cpp'), ...
        fullfile(gpudir,'CudaGPU.cpp'), ...
        fullfile(gpudir,'AbstractInterface.cpp'), ...
        fullfile(gpudir,'MatlabInterface.cpp'), ...
        fullfile(gpudir,'Lattice.cpp'), ...
        fullfile(gpudir,'PassMethodFactory.cpp'), ...
        fullfile(gpudir,'SymplecticIntegrator.cpp'), ...
        fullfile(gpudir,'IdentityPass.cpp'), ...
        fullfile(gpudir,'DriftPass.cpp'), ...
        fullfile(gpudir,'StrMPoleSymplectic4Pass.cpp'), ...
        fullfile(gpudir,'BndMPoleSymplectic4Pass.cpp'), ...
        fullfile(gpudir,'StrMPoleSymplectic4RadPass.cpp'), ...
        fullfile(gpudir,'BndMPoleSymplectic4RadPass.cpp'), ...
        fullfile(gpudir,'CavityPass.cpp'), ...
        fullfile(gpudir,'RFCavityPass.cpp'), ...
        fullfile(gpudir,'ExactDriftPass.cpp'), ...
        fullfile(gpudir,'ExactMultipolePass.cpp'), ...
        fullfile(gpudir,'ExactMultipoleRadPass.cpp'), ...
        fullfile(gpudir,'ThinMPolePass.cpp'), ...
        '-lcuda','-lnvrtc');
end

[warnmess,warnid]=lastwarn; %#ok<ASGLU>
if strcmp(warnid,'MATLAB:mex:GccVersion_link')
    warning('Disabling the compiler warning');
end
oldwarns=warning('OFF','MATLAB:mex:GccVersion_link');

% Diffusion matrices
cdir=fullfile(atroot,'atphysics','Radiation');
compile([alloptions, {passinclude}], fullfile(cdir,'findmpoleraddiffmatrix.c'));

% RDTs
cdir=fullfile(atroot,'atphysics','NonLinearDynamics');
compile(alloptions, fullfile(cdir,'RDTelegantAT.cpp'));

% NAFF
cdir=fullfile(atroot,'atphysics','nafflib');
compile(alloptions, fullfile(cdir,'nafflib.c'),...
                    fullfile(cdir,'modnaff.c'),...
                    fullfile(cdir,'complexe.c'));

% Find all files matching '*Pass.c' wildcard
cfiles = dir(fullfile(pdir,'*Pass.c'));
passmethods = {cfiles.name};
if ~c_only
    % Find all files matching '*Pass.cc' wildcard
    ccfiles = dir(fullfile(pdir, '*Pass.cc'));
    passmethods = [passmethods ccfiles.name];
end
% Eliminate invisible files
ok=cellfun(@(nm) nm(1)~='.',passmethods,'UniformOutput',false);
passmethods = passmethods(cell2mat(ok));
try
    generate_passlist(pdir,passmethods);
catch err
    fprintf(2,'\nCannot generate the list of passmethods: %s\n\n', err.message);
end

for i = 1:length(passmethods)
    PM = fullfile(pdir,[passmethods{i}]);
    if exist(PM,'file')
        try
            if exist('map2','var')
                try
                    compile([map1, alloptions], PM);
                catch
                    compile([map2, alloptions], PM);
                end
            else
                compile([map1, alloptions], PM);
            end
        catch err
            if fail
                rethrow(err);
            else
                fprintf(2, 'Could not compile %s:\n%s\n', PM, err.message);
            end
        end
    else
        if fail
            error('AT:atmexall', '%s not found\n', PM);
        else
            fprintf(2,'%s not found, skip to next\n', PM);
        end
    end
end

warning(oldwarns.state,oldwarns.identifier);

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

	function flags=pc_flags(flags, cflags,ldflags)
        if ~isempty(ldflags)
            flags=[{sprintf('LINKFLAGS="$LINKFLAGS %s"',strjoin(ldflags))}, flags];
        end
        if ~isempty(cflags)
            flags=[{sprintf('COMPFLAGS="$COMPFLAGS %s"',strjoin(cflags))}, flags];
        end
    end
        

    function flags=unix_flags(flags, cflags, ldflags, exportarg, map)
        if nargin > 3
            mapfile=fullfile(pdir, map);
            if verLessThan('matlab','8.3')                          %  R2008a < < R2014a
                exp=sprintf([exportarg '"%s"'], mapfile);
                def_ldflags=regexprep(mex.getCompilerConfigurations('C').Details.LinkerFlags,['\s(' exportarg ')([^\s,]+)'],'');
                flags=[{sprintf('LDFLAGS=%s',strjoin([{def_ldflags}, ldflags, {exp}]))}, flags];
                ldflags = {};
            elseif verLessThan('matlab','9.1')                      %           < R2016b
                flags=[{sprintf('LINKEXPORT=%s"%s"',exportarg,mapfile)}, flags];
            else                                                    %           >= R2016b
                %           Starting from R2016b, Matlab introduced a new entry point in MEX-files
                %           The "*.mapext" files defines this new entry point
                flags=[{sprintf('VERSIONMAP="''%sext''"',mapfile),...
                    sprintf('FUNCTIONMAP="''%s''"',mapfile)}, flags];
            end
        end
        if ~isempty(ldflags)
            flags = [{sprintf('LDFLAGS="$LDFLAGS %s"',strjoin(ldflags))} flags];
        end
        if ~isempty(cflags)
            flags=[{sprintf('CFLAGS="$CFLAGS %s"',strjoin(cflags))}, flags];
        end
    end

    function compile(mexargs, varargin)
        [fpath, fname, ~] = fileparts(varargin{1});
        target = strjoin({fullfile(fpath, fname), mexext}, '.');
        if force || ~exist(target, 'file') || ...
                any(cellfun(@(f) getdate(f) >  getdate(target), varargin))
%               any(cellfun(@(f) dir(f).datenum >  dir(target).datenum, varargin))  % Not accepted in R2016b
            disp(['mex ',strjoin([mexargs, {'-outdir', fpath}, varargin])]);
            mex(mexargs{:},'-outdir', fpath, varargin{:});
        end
        
        function d=getdate(f)
            dt=dir(f);
            d=dt.datenum;
        end
    end

    function include=select_omp()
        if exist('/usr/local/include/omp.h', 'file')            % Homebrew
            include='-I/usr/local/include';
        elseif exist('/opt/local/include/libomp/omp.h', 'file') % MacPorts
            include='-I/opt/local/include/libomp';
        else
            error('AT:MissingLibrary', strjoin({'', ...
                'libomp.dylib must be installed with your favourite package manager:', '', ...
                'Use "$ brew install libomp"', ...
                'Or  "$ sudo port install libomp"'...
                }, '\n'));
        end
    end

    % Generate "element.gpuh" from "Element.h" fpr GPU extensions
    function generate_gpuheader(gpudir)

        disp(strjoin(['Generating ',fullfile(gpudir,"element.gpuh")]));
        lines = [];
        f = fopen(fullfile(gpudir,"Element.h"),'r');
        line = fgetl(f);
        while ischar(line)
            lines = [lines;convertCharsToStrings(line)];
            line = fgetl(f);
        end
        fclose(f);

        % Add type definitions
        idx = startsWith(lines,'// DEFINE_TYPE');
        I = find(idx>0,1);
        if numel(I)==0
            error('Invalid Element.h file "// DEFINE_TYPE" not found');
        end
        types = ["typedef signed char        int8_t;";...
                 "typedef signed short       int16_t;";...
                 "typedef signed int         int32_t;";...
                 "typedef signed long long   int64_t;";...
                 "typedef unsigned char      uint8_t;";...
                 "typedef unsigned short     uint16_t;";...
                 "typedef unsigned int       uint32_t;";...
                 "typedef unsigned long long uint64_t;"];
        lines = [ lines(1:I-1);types;lines(I+1:end) ];

        %insert raw string marker
        lines = ["R""(";lines];
        lines = [lines;")"""];

        % Write the file
        f = fopen(fullfile(gpudir,"element.gpuh"),'w');
        fprintf(f,"%s\n",lines);
        fclose(f);

    end

end
