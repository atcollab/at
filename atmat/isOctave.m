function retval = isOctave
%isOctave check if running Octave
%
%  [RETVAL]=ISOCTAVE()
%    Check if running Octave
%
%  OUTPUTS
%  1. RETVAL          boolean is running Octave
%
  retval = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
