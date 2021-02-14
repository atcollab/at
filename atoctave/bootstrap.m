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
  INCLUDES=strjoin({'-I../atintegrators'}, ' ');
  CFLAGS=strjoin({'-DOCTAVE', ...
                  '-DMATLAB_MEX_FILE', ...
                  '-Wall', ...
                  '-Wextra', ...
                  '-Wno-unused-function', ...
                  '-Wno-unused-parameter'}, ' ');
  if debugMode
    CFLAGS=strjoin({CFLAGS, ...
                    '-ggdb3', ...
                    '-O0'}, ' ');
  else
    CFLAGS=strjoin({CFLAGS, ...
                    '-O3'}, ' ');
  endif
  if useOpenMP
    if ismac
      OMPFLAGS='-Xpreprocessor -fopenmp -lomp';
    else
      OMPFLAGS='-fopenmp';
    endif
  else
    OMPFLAGS='';
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
  %% get source files

  %% search for all *Pass.c and *Pass.cpp in atintegrators folder
  atintegrators = glob(fullfile(atroot, '..', 'atintegrators', '*Pass.c*'));
  for i = 1:length(atintegrators)
    source = atintegrators{i};
    target = replaceext(source, '.mex');
    atcompilemex(target, strjoin({CFLAGS, OMPFLAGS}), INCLUDES, {source});
  endfor

  nonlineardynamics = glob(fullfile(atroot, 'atphysics', 'NonLinearDynamics', '*.c*'));
  atcompilemex(replaceext(nonlineardynamics{1}, '.mex'), CFLAGS, "", nonlineardynamics);

  radiation = glob(fullfile(atroot, 'atphysics', 'Radiation', '*.c*'));
  atcompilemex(replaceext(radiation{1}, '.mex'), CFLAGS, INCLUDES, radiation);

  attrack = glob(fullfile(atroot, 'attrack', '*.c'));
  atcompilemex(replaceext(attrack{1}, '.mex'), strjoin({CFLAGS, OMPFLAGS}), INCLUDES, attrack);

  nafflib = glob(fullfile(atroot, 'atphysics', 'nafflib', '*.c'));
  nafflib = nafflib(~strcmp(nafflib, fullfile(atroot, 'atphysics', 'nafflib', 'example.c')));
  atcompilemex(fullfile(atroot, 'atphysics', 'nafflib', 'nafflib.mex'), CFLAGS, INCLUDES, nafflib);

  %% }}}

  %% {{{ Functions
  function atcompilemex(target, CFLAGS, INCLUDES,sources)
    %% if file does not exist or is older then sources
    if ~exist(target, 'file') || ...
      any(arrayfun(@(f) stat(f{1}).mtime, sources) > stat(target).mtime)
      mexcommand = strjoin({'mex', ...
                            CFLAGS, ...
                            INCLUDES, ...
                            '-o', ...
                            target, ...
                            strjoin(sources, ' ')},
                           ' ');
      disp(mexcommand);
      eval(mexcommand);
    endif
  end

  function output=replaceext(string, ext)
    [~, ~, ~, m] = regexp(string, "(?i)\\.[a-z0-9]+$");
    if length(m) > 0
      output = strrep(string, m{1}, ext);
    else
      output = string;
    end
  end
  %% }}}
end
