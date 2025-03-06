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
                    rethrow(err);
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
end
