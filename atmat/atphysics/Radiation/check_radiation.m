function varargout=check_radiation(ring,varargin)
%CHECK_RADIATION	Obsolete: check the radiation state of a ring
%
% Kept for compatibility> The function name is misleading, because the
% function checks not only the presence of synchrotron radiation, but more
% generally of all elements modifying the longitudinal momentum.
%
% <a href="matlab:help check_6d">check_6d</a> is an exact copy of this function and should preferably be
% used.
%
% See also: CHECK_6D, ATGETRINGPROPERTIES, ATENABLE_6D, ATDISABLE_6D

[varargout{1:nargout}]=check_6d(ring,varargin{:});

end
