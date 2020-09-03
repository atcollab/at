function varargout=GWigSymplecticRadPass(varargin) %#ok<STOUT>
%GWigSymplecticRadPass - generic symplectic integrator for wigglers
%
% Similar to GWigSymplecticPass but with radiation included
%
% For documentation on the integrator, see GWigSymplecticPass
%
%see also: GWigSymplecticRadPass

error('AT:MissingMexFile','"%s.%s" is missing. Did you run "atmexall" ?',...
    mfilename,mexext);
end
