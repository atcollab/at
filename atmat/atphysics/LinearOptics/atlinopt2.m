function varargout = atlinopt2(ring,varargin)
%ATLINOPT2 Performs the linear analysis of UNCOUPLED lattices
%
%[RINGDATA,ELEMDATA] = ATLINOPT2(RING,REFPTS)
%
% IMPORTANT!!! ATLINOPT2 assumes a constant momentum deviation.
%   PassMethods used for any element in the RING SHOULD NOT
%   1.change the longitudinal momentum dP
%     (cavities , magnets with radiation, ...)
%   2.have any time dependence (localized impedance, fast kickers, ...)
%
%RINGDATA is a structure array with fields:
%   tune          1x2 tunes
%   chromaticity  1x2 chromaticities (only with get_chrom or get_w flags)
%
%ELEMDATA is a structure array with fields:
%   SPos        - longitudinal position [m]
%   ClosedOrbit - 4x1 closed orbit vector with components
%                 x, px, y, py (momentums, NOT angles)
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   M           - 4x4 transfer matrix M from the beginning of RING
%                 to the entrance of the element [2]
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%   mu          - [mux, muy] horizontal and vertical betatron phase advances
%   W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
%                           get_w flag)
%
%   All values are specified at the entrance of each element specified in REFPTS.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
%   See further explanation of REFPTS in the 'help' for FINDSPOS
% 
%   Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
%   Example: 
%   >> [ringdata, elemdata] = ATLINOPT2(ring,1:length(ring));
%   >> beta = cat(1,elemdata.beta);
%   >> s = cat(1,elemdata.SPos);
%   >> plot(S,beta)
%
% [...] = ATLINOPT2(...,'get_chrom')
%   Trigger the computation of chromaticities
%
% [...] = ATLINOPT2(...,'get_w')
%   Trigger the computation of chromatic amplitude functions (time consuming)
%
% [...] = ATLINOPT2(...,'dp',DPP)
%   Analyses the off-momentum lattice by specifying the central
%   off-momentum value
%
% [...] = ATLINOPT2(...,'ct',DCT)
%   Analyses the off-momentum lattice by specifying the path lengthening
%   corresponding to a frequency shift. The resulting deltap/p is returned
%   in the 5th component of the ClosedOrbit field.
%
% [...] = ATLINOPT2(...,'orbit',ORBITIN)
%   Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
%   of initial conditions is used: [x0; px0; y0; py0; DP; 0].
%   The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit if RING is not a
%   ring or to avoid recomputing the closed orbit if is already known.
%
% [...] = ATLINOPT2(...,'twiss_in',TWISSIN)
%   Computes the optics for a transfer line.
%
% TWISSIN is a scalar structure with fields:
%   ClosedOrbit - 4x1 initial closed orbit. Default: zeros(4,1)
%   Dispersion  - 4x1 initial dispersion.   Default: zeros(4,1)
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax0, betay0] vector
%   alpha       - [alphax0, alphay0] vector
%
%  REFERENCES
%	[1] D.Edwards,L.Teng IEEE Trans.Nucl.Sci. NS-20, No.3, p.885-888, 1973
%	[2] E.Courant, H.Snyder
%	[3] Brian W. Montague Report LEP Note 165, CERN, 1979
%
%  See also atlinopt4, atlinopt6

[varargout{1:nargout}]=atlinopt4(ring,varargin{:},'coupled',false);
end
