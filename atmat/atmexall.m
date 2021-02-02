function atmexall(varargin)
%ATMEXALL builds all AT platform dependent mex-files from C-sources
%
% On UNIX platform, the GNU gcc compiler must be installed and
% properly configured.
% On Windows, Microsoft Visual C++ is required

if isOctave
  error('atmexall does not work in Octave, use atoctave/bootstrap instead')
end

[~,varargs]=getflag(varargin,'-openmp');

atoptions={['-D',computer]};

if isunix && ~ismac
    LIBDL={'-ldl'};
else
    LIBDL={};
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

passinclude = ['-I' fullfile(atroot,'..','atintegrators','')];
alloptions=[atoptions varargs];

% atpass
cdir=fullfile(atroot,'attrack','');
mexargs=[alloptions, {passinclude}, LIBDL, {'-outdir', cdir, fullfile(cdir,'atpass.c')}];
disp(['mex ',strjoin(mexargs)]);
mex(mexargs{:});

[warnmess,warnid]=lastwarn; %#ok<ASGLU>
if strcmp(warnid,'MATLAB:mex:GccVersion_link')
    warning('Disabling the compiler warning');
end
oldwarns=warning('OFF','MATLAB:mex:GccVersion_link');

% Diffusion matrices
cdir=fullfile(atroot,'atphysics','Radiation');
mexargs=[alloptions,{passinclude, '-outdir', cdir, fullfile(cdir,'findmpoleraddiffmatrix.c')}];
disp(['mex ',strjoin(mexargs)]);
mex(mexargs{:});

% RDTs
cdir=fullfile(atroot,'atphysics','NonLinearDynamics');
mexargs=[alloptions,{'-outdir', cdir, fullfile(cdir,'RDTelegantAT.cpp')}];
disp(['mex ',strjoin(mexargs)]);
mex(mexargs{:});

% NAFF
cdir=fullfile(atroot,'atphysics','nafflib');
mexargs=[alloptions,{'-outdir', cdir, fullfile(cdir,'nafflib.c'),...
                                      fullfile(cdir,'modnaff.c'),...
                                      fullfile(cdir,'complexe.c')}];
disp(['mex ',strjoin(mexargs)]);
mex(mexargs{:});


% Navigate to the directory that contains pass-methods 
mexpassmethod('all',atoptions{:},varargin{:});
warning(oldwarns.state,oldwarns.identifier);

% ADD 'MEXING' instructions for other C files
%disp('ALL mex-files created successfully')
end
