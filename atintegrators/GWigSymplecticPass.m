function varargout=GWigSymplecticPass(varargin) %#ok<STOUT>
%GWigSymplecticPass - generic symplectic integrator for wigglers
% 
% Integrator based on:
%
%   "Explicit symplectic integrator for s-dependent static magnetic field"
%   by Wu, Y. K. and Forest, E. and Robin, D. S.
%   PhysRevE.68.046502
%   https://link.aps.org/doi/10.1103/PhysRevE.68.046502
%
%The field of an horizontal wiggler is represented by a sum of harmonics,
%each being described by a 6x1 column vector B such as:
%
%   Bx/Bmax =  B2 * B3/B4 * sin(B3*kw*x) sinh(B4*kw*z) cos(B5*kw*s + B6)
%   Bz/Bmax = -B2 * cos(B3*kw*x) cosh(B4*kw*z) cos(B5*kw*s + B6)
%   Bs/Bmax =  B2 * B5/B4 * cos(B3*kw*x) sinh(B4*kw*z) sin(B5*kw*s + B6)
%
%   with kw = 2*Pi/Lw   and   B4^2 = B3^2 + B5^2
%
% Bmax: Peak wiggler field
% Lw:   Period length
% B2:   Fraction of Bmax for the given harmonic
% B3:   Field dependence on x
% B4:   Field dependence on z
% B5:   Harmonic number
% B6:   Longitudinal phase of the given harmonic
%
%So the horizontal wiggler component is described by a 6 x NHharm matrix 'By'
%
%As an example, the fundamental of an infinitely wide horizontal wiggler
%(no horizontal dependence of the vertical field) is represented by:
%
%                By=[1;1;0;1;1;0]
%
%Similarly, the field of an harmonic of a vertical wiggler is:
%
%   Bx/Bmax =  B2 * cosh(B3*kw*x) cos(B4*kw*z) cos(B5*kw*s + B6)
%   Bz/Bmax = -B2 * B4/B3 * sinh(B3*kw*x) sin(B4*kw*y) cos(B5*kw*s + B6)
%   Bs/Bmax = -B2 * B5/B3 * sinh(B3*kw*x) cos(B4*kw*z) sin(B5*kw*s + B6)
%
%   with kw = 2*Pi/Lw   and   B3^2 = B4^2 + B5^2
%
%And the vertical wiggler component is described by a 6 x NVharm matrix 'Bx'
%
%An arbitralily polarized wiggler is described be the combination of
%horizontal and vertical wiggler fields with an adjustable phasing.
%
%see also: GWigSymplecticRadPass

error('AT:MissingMexFile','"%s.%s" is missing. Did you run "atmexall" ?',...
    mfilename,mexext);
end
