function atmexall(varargin)
%ATMEXALL build all AT platform dependent mex-files from C-sources
%
%ATMEXALL -openmp       Build the integrators for OpenMP parallelisation
%
%ATMEXALL -only_new     Build only the outdated components
%
%ATMEXALL [options]     Other options are transmitted to the mex command
%
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required

pdir=fullfile(fileparts(atroot),'atintegrators');
[openmp,varargs]=getflag(varargin,'-openmp');
[miss_only,varargs]=getflag(varargs,'-missing');
[debug,varargs]=getflag(varargs,'-g');
force=~miss_only;

cpt = computer('arch');
if ispc()
    computer = 'PCWIN64';
elseif strncmp(cpt, 'darwin', 6)
    computer = 'MACI64';
elseif strncmp(cpt, 'gnu-linux)', 9)
    computer = 'GLNXA64';
else
    computer = 'UNKNOWN';
end

atoptions={'-DOCTAVE', ...
    '-DMATLAB_MEX_FILE', ...
    ['-D' computer]};

if debug
    atoptions=[atoptions {'-ggdb3', '-O0'}];
else
    atoptions=[atoptions {'-O3'}];
end

if openmp
    if ispc()
        ompoptions={'/openmp'};
    elseif ismac()
        ompoptions={'-I/usr/local/include','-Xpreprocessor','-fopenmp','-lomp'};
    else
        ompoptions={'-fopenmp','-lgomp'};
    end
else
    ompoptions={};
end

passinclude = ['-I' pdir];
alloptions=[atoptions varargs];

% atpass
cdir=fullfile(atroot,'attrack','');
compile([alloptions, {passinclude}, ompoptions], fullfile(cdir,'atpass.c'));

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
            compile([alloptions, ompoptions], PM);
        catch err
            fprintf(2, 'Could not compile %s: %s\n', PM, err.message);
        end
    else
        fprintf(2,'%s not found, skip to next\n', PM);
    end
end

    function compile(mexargs, varargin)
        [fpath, fname, ~] = fileparts(varargin{1});
        target = strjoin({fullfile(fpath, fname), mexext}, '.');
        if force || ~exist(target, 'file') || ...
                any(cellfun(@(f) stat(f).mtime > stat(target).mtime, varargin));
            disp(['mex ',strjoin([mexargs, {'-o', target}, varargin])]);
            mex(mexargs{:},'-o', target, varargin{:});
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
