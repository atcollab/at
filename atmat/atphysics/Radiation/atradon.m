function varargout=atradon(ring, varargin)
%ATRADON    Obsolete: switches RF and radiation on
%
% Kept for compatibility. The function name is misleading, because the
% function acts not only on synchrotron radiation, but more generally on
% all elements modifying the longitudinal momentum.
%
% <a href="matlab:help atenable_6d">atenable_6d</a> is an exact copy of this function and should preferably be
% used.
%
%  See also ATENABLE_6D, ATDISABLE_6D, CHECK_6D, ATCAVITYON, ATCAVITYOFF

[varargout{1:nargout}]=atenable_6d(ring, varargin{:});

end
