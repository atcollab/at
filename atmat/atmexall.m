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
%	-opencl OCL_PATH Build the GPU tracking support using OpenCL
%               Use "-opencl default" for using standard OpenCL install
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
[opencl,varargs]=getoption(varargs,'-opencl','None');
[miss_only,varargs]=getflag(varargs,'-missing');
[c_only,varargs]=getflag(varargs,'-c_only');
[fail,varargs]=getflag(varargs,'-fail');
force=~miss_only;

atoptions=[varargs {['-D',computer]}];

if isunix && ~ismac
    LIBDL={'-ldl'};
else
    LIBDL={};
end

if isunix
    cflags={'-std=c99'};
    compformat='CFLAGS="$CFLAGS %s"';
    linkformat='LDFLAGS="$LDFLAGS %s"';
else
    cflags={};
    compformat='COMPFLAGS="$COMPFLAGS %s"';
    linkformat='LINKFLAGS="$LINKFLAGS %s"';
end
cppflags={};
ldflags={};

if openmp
    if ispc()
        cflags=[cflags {'/openmp'}];
        ompflags={};
    elseif ismac()
        cflags=[cflags {'-Xpreprocessor', '-fopenmp'}];
        ompflags = select_omp();
    else
        cflags=[cflags {'-fopenmp'}];
        ompflags={...
            sprintf('-L"%s"',fullfile(matlabroot,'sys','os',computer('arch'))),...
            '-liomp5'};
    end
else
    ompflags={};
end

mexflags = make_flags(cflags, ldflags);
mexcpp = make_flags(cppflags, ldflags);
if ispc()
    mapc=mexflags;
    mapcpp=mexcpp;
elseif ismac()
    mapc=make_flags(cflags, [ldflags {'-Wl,-exported_symbol,_trackFunction'}]);
    mapcpp=make_flags(cppflags, [ldflags {'-Wl,-exported_symbol,_trackFunction'}]);
else
    exportflag = sprintf('VERSIONMAP="%s"',fullfile(pdir, 'mexFunctionGLNX86.mapext'));
    mapc=[mexflags {exportflag, 'LINKEXPORT=""'}];
    mapcpp=[mexcpp {exportflag, 'LINKEXPORT=""'}];
end

if ~verLessThan('matlab','9.4')         % >= R2018a
    atoptions = [atoptions,{'-R2018a'}];
elseif ~verLessThan('matlab','7.11')	% >= R2010b
    atoptions = [atoptions,{'-largeArrayDims'}];
end
if ~verLessThan('matlab','8.3')         % >= R2014a
    atoptions = [atoptions,{'-silent'}];
end

lastwarn('');
passinclude = ['-I' pdir];

% atpass
cdir=fullfile(atroot,'attrack','');
compile([atoptions, {passinclude}, LIBDL, mexflags, ompflags], fullfile(cdir,'atpass.c'));
compile([atoptions, mexflags, ompflags],fullfile(cdir,'coptions.c'))

% gpuextensions
iscuda = ~strcmp(cuda,'None');
isopencl = ~strcmp(opencl,'None');
if iscuda || isopencl
    gpudir=fullfile(fileparts(atroot),'atgpu','');
    generate_gpuheader(gpudir)
    if ispc()
        % TODO
        error('AT:atmexall', 'GPU windows not supported');
    elseif ismac()
        % TODO
        error('AT:atmexall', 'GPU ismac not supported');
    else
        if iscuda
            gpuflags = {sprintf('-I"%s"',gpudir),...
                        sprintf('-I"%s/include"',cuda),...
                        sprintf('-L"%s/lib64"',cuda),...
                        sprintf('LDFLAGS=$LDFLAGS -Wl,-rpath,"%s/lib64"',cuda),...
                        '-DCUDA'};
        else
            if strcmp(opencl,'default')
                gpuflags = {sprintf('-I"%s"',gpudir),...
                            '-DOPENCL'};
            else
                gpuflags = {sprintf('-I"%s"',gpudir),...
                            sprintf('-I"%s/include"',opencl),...
                            sprintf('-L"%s/lib"',opencl),...
                            sprintf('LDFLAGS=$LDFLAGS -Wl,-rpath,"%s/lib64"',opencl),...
                            '-DOPENCL'};
            end
        end
    end

    % gpuinfo
    if iscuda
        compile([alloptions, {passinclude}, gpuflags], ...
            fullfile(cdir,'gpuinfo.cpp'),...
            fullfile(gpudir,'MatlabInterface.cpp'), ...
            fullfile(gpudir,'AbstractInterface.cpp'), ...
            fullfile(gpudir,'CudaGPU.cpp'), ...
            fullfile(gpudir,'AbstractGPU.cpp'), ...
            '-lcuda','-lnvrtc');
    else
        compile([alloptions, {passinclude}, gpuflags], ...
            fullfile(cdir,'gpuinfo.cpp'),...
            fullfile(gpudir,'MatlabInterface.cpp'), ...
            fullfile(gpudir,'AbstractInterface.cpp'), ...
            fullfile(gpudir,'OpenCLGPU.cpp'), ...
            fullfile(gpudir,'AbstractGPU.cpp'), ...
            '-lOpenCL');
    end

    % gpupass
    if iscuda
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
            fullfile(gpudir,'CorrectorPass.cpp'), ...
            fullfile(gpudir,'AperturePass.cpp'), ...
            '-lcuda','-lnvrtc');
    else
        compile([alloptions, {passinclude}, gpuflags], ...
            fullfile(cdir,'gpupass.cpp'),...
            fullfile(gpudir,'AbstractGPU.cpp'), ...
            fullfile(gpudir,'OpenCLGPU.cpp'), ...
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
            fullfile(gpudir,'CorrectorPass.cpp'), ...
            fullfile(gpudir,'AperturePass.cpp'), ...
            '-lOpenCL');                 
    end

end

[warnmess,warnid]=lastwarn; %#ok<ASGLU>
if strcmp(warnid,'MATLAB:mex:GccVersion_link')
    warning('Disabling the compiler warning');
end
oldwarns=warning('OFF','MATLAB:mex:GccVersion_link');

% Diffusion matrices
cdir=fullfile(atroot,'atphysics','Radiation');
compile([atoptions, {passinclude}, mexflags], fullfile(cdir,'findmpoleraddiffmatrix.c'));
compile([atoptions, {passinclude}, mexflags], fullfile(cdir,'FDW.c'));
compile([atoptions, {passinclude}, LIBDL, mexflags], fullfile(cdir,'diffusion_matrix.c'));

% RDTs
cdir=fullfile(atroot,'atphysics','NonLinearDynamics');
compile([atoptions, mexcpp], fullfile(cdir,'RDTelegantAT.cpp'));

% NAFF
cdir=fullfile(atroot,'atphysics','nafflib');
compile([atoptions, mexflags], fullfile(cdir,'nafflib.c'),...
                               fullfile(cdir,'modnaff.c'),...
                               fullfile(cdir,'complexe.c'));

% Passmethods
cfiles = passcomp(pdir, '*Pass.c',[atoptions, mapc, ompflags]);
if ~c_only
    % Compile '*Pass.cc' methods
    ccfiles = passcomp(pdir, '*Pass.cc',[atoptions, mapcpp, ompflags]);
else
    ccfiles={};
end

try
    generate_passlist(pdir,[cfiles ccfiles]);
catch err
    fprintf(2,'\nCannot generate the list of passmethods: %s\n\n', err.message);
end

warning(oldwarns.state,oldwarns.identifier);

    function meths=passcomp(pdir, pattern, compargs)
        meths = dir(fullfile(pdir,pattern));
        % Keep the names
        meths = {meths.name};
        % Discard hidden files
        meths = meths(cellfun(@(nm) nm(1)~='.',meths));
        for im=1:length(meths)
            pmeth = fullfile(pdir,[meths{im}]);
            try
                compile(compargs, pmeth);
            catch errcomp
                if fail
                    rethrow(errcomp);
                else
                    fprintf(2, 'Could not compile %s:\n%s\n', pmeth, errcomp.message);
                end
            end
        end
        % Remove extension
        meths = cellfun(@dpl,meths, 'UniformOutput',false);

        function nm=dpl(fname)
            [~,nm,~] = fileparts(fname);
        end
    end


    function generate_passlist(pdir,passmethodnames)
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

    function flags=make_flags(cflags, ldflags)
        flags = {};
        if ~isempty(ldflags)
            flags = [{sprintf(linkformat,strjoin(ldflags))} flags];
        end
        if ~isempty(cflags)
            flags = [{sprintf(compformat,strjoin(cflags))} flags];
        end
    end

    function compile(mexargs, varargin)
        [fpath, fname, ~] = fileparts(varargin{1});
        target = strjoin({fullfile(fpath, fname), mexext}, '.');
        if force || ~exist(target, 'file') || ...
                any(cellfun(@(f) getdate(f) >  getdate(target), varargin))
            disp(['mex ',strjoin([mexargs, {'-outdir', fpath}, varargin])]);
            mex(mexargs{:},'-outdir', fpath, varargin{:});
        end
        
        function d=getdate(f)
            dt=dir(f);
            d=dt.datenum;
        end
    end

    function flags=select_omp()
        arch=computer('arch');
        if strcmp(arch, 'maca64')
            homeb='/opt/homebrew/opt/libomp/include';
            libdir=fullfile(matlabroot,'bin',arch);
            libname='omp';
        else
            homeb='/usr/local/include';
            libdir=fullfile(matlabroot,'sys','os',arch);
            libname='iomp5';
        end
        if exist([homeb '/omp.h'], 'file')                      % Homebrew
            incdir=homeb;
        elseif exist('/opt/local/include/libomp/omp.h', 'file') % MacPorts
            incdir='/opt/local/include/libomp';
        else
            error('AT:MissingLibrary', strjoin({'', ...
                'libomp.dylib must be installed with your favourite package manager:', '', ...
                'Use "$ brew install libomp"', ...
                'Or  "$ sudo port install libomp"'...
                }, '\n'));
        end
        flags={['-I' incdir], sprintf('-L"%s"',libdir),['-l' libname]};
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
