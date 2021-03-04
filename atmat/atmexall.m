function atmexall(varargin)
%ATMEXALL build all AT platform dependent mex-files from C-sources
%
%ATMEXALL option1 ... optionN
%
% AT Options:
%
%	-missing	Build only the outdated components
%	-openmp     Build the integrators for OpenMP parallelisation
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
        cflags={'-Xpreprocessor', '-fopenmp'};
        ompflags={'-I/usr/local/include',...
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
            fprintf(2, 'Could not compile %s: %s\n', PM, err.message);
        end
    else
        fprintf(2,'%s not found, skip to next\n', PM);
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
                flags=[{sprintf('VERSIONMAP="%sext"',mapfile),...
                    sprintf('FUNCTIONMAP="%s"',mapfile)}, flags];
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
                any(cellfun(@(f) dir(f).datenum >  dir(target).datenum, varargin))
            disp(['mex ',strjoin([mexargs, {'-outdir', fpath}, varargin])]);
            mex(mexargs{:},'-outdir', fpath, varargin{:});
        end
    end
end
