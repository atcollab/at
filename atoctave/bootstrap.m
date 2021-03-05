function bootstrap(debugMode=false,useOpenMP=false)
%bootstrap prepare Octave for AT
%
%  ISOCTAVE(DEBUGMODE)
%    Prepare Octave for AT
%
%  INPUTS
%  1. DEBUGMODE          build debug versions of mex files. defaults to false
%  1. USEOPENMP          parallelize C loops in intergrators. defaults to true
%

  %% {{{ Constants

  %% for debugging use flags -O0 -ggdb3
  CFLAGS={'-missing', '-Wall', '-Wextra', ...
          '-Wno-unused-function', ...
          '-Wno-unused-parameter'};
  if debugMode
    CFLAGS=[CFLAGS {'-g'}];
  endif
  if useOpenMP
    CFLAGS=[CFLAGS {'-openmp'}];
  endif

  %% }}}

  %% {{{ Prepare paths
  [octavedir, ~, ~] = fileparts(which(mfilename));
  run(fullfile(octavedir,'..','atmat','atpath.m'));
  atroot = atroot();
  disp(octavedir);
  addpath(genpath(octavedir));
  %% }}}

  %% {{{ Check Octave
  if ~isOctave
    error("This script is for Octave only");
  end

  if octaveVersion < 6
    error("Minimal supported Octave version is 6. See https://savannah.gnu.org/bugs/?39257");
  end
  %% }}}

  %% {{{ Pkg config
  %% TODO: check that these paths work with windows
  try
    pkg load optim;
  catch
    %% install from source forge
    display('Installing required packages...');
    pkg('prefix', ...
        fullfile('~', '.octave'), ...
        fullfile('~', '.octave_arch'));
    pkg('install', '-forge', '-local', 'struct');
    pkg('install', '-forge', '-local', 'io');
    pkg('install', '-forge', '-local', 'statistics');
    pkg('install', '-forge', '-local', 'optim');
    display('Done');
    pkg load optim;
  end_try_catch
  %% }}}

  %% {{{ Compile MEX files
  atmexall(CFLAGS{:});
  %% }}}
end
