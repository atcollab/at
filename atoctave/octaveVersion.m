function [maj, min] = octaveVersion
%octaveVersion get major and minor version numbers
%
%  [MAJ,MIN]=OCTAVEVERSION()
%    Get major and minor version numbers
%
%  OUTPUTS
%  1. MAJ          major version number
%  2. MIN          minor version number
%
  retval = strsplit(OCTAVE_VERSION, '.');
  maj = str2num(retval{1});
  min = str2num(retval{2});
end
