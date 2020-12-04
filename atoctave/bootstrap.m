function bootstrap(debug=false)
%bootstrap prepare Octave for AT
%
%  ISOCTAVE(DEBUG)
%    Prepare Octave for AT
%
%  INPUTS
%  1. DEBUG          build debug versions of mex files. defaults to false
%

  %% {{{ Constants

  %% for debugging use flags -O0 -ggdb3
  INCLUDES=strjoin({'-I../atintegrators'}, ' ');
  CFLAGS=strjoin({'-DOCTAVE', ...
                  '-DAT_MODE', ...
                  '-DMATLAB_MEX_FILE', ...
                  '-Wall', ...
                  '-Wextra', ...
                  '-Wno-unused-function', ...
                  '-Wno-unused-parameter'}, ' ');
  if debug
    CFLAGS=strjoin({CFLAGS, ...
                    '-ggdb3', ...
                    '-O0'}, ' ');
  else
    CFLAGS=strjoin({CFLAGS, ...
                    '-O3'}, ' ');
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
  atintegrators = glob(fullfile(atroot, '..', 'atintegrators', '*Pass.c*'));
  nonlineardynamics = glob(fullfile(atroot, 'atphysics', 'NonLinearDynamics', '*.c'));
  radiation = glob(fullfile(atroot, 'atphysics', 'Radiation', '*.c'));
  attrack = glob(fullfile(atroot, 'attrack', '*.c'));
  nafflib = glob(fullfile(atroot, 'atphysics', 'nafflib', '*.c'));
  nafflib = nafflib(~strcmp(nafflib, fullfile(atroot, 'atphysics', 'nafflib', 'example.c')));

  sources = [atintegrators; nonlineardynamics; nafflib; radiation; attrack];

  for i = 1:length(sources)
    source = sources{i};
    if endsWith(source, '.cc')
      target = strrep(source, '.cc', '.mex');
    else
      target = strrep(source, '.c', '.mex');
    endif
    if ~exist(target, 'file')
      mexcommand = strjoin({'mex', ...
                            CFLAGS, ...
                            INCLUDES, ...
                            '-o', ...
                            target, ...
                            source}, ' ');
      disp(mexcommand);
      eval(mexcommand);
    endif
  endfor
  %% }}}

end
