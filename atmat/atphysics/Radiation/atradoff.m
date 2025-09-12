function varargout=atradoff(ring, varargin)
%ATRADOFF   Obsolete: switches RF and radiation off
%
% Kept for compatibility. The function name is misleading, because the
% function acts not only on synchrotron radiation, but more generally on
% all elements modifying the longitudinal momentum.
%
% <a href="matlab:help atdisable_6d">atdisable_6d</a> is an exact copy of this function and should preferably be
% used.
%
%  See also ATDISABLE_6D, ATENABLE_6D, CHECK_6D, ATCAVITYOFF, ATCAVITYON

[varargout{1:nargout}]=atdisable_6d(ring, varargin{:});
    
end
