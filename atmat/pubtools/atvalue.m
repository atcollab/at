function v=atvalue(lindata,code,idx)
%ATVALUE			extract array from lindata structure
%
%V=ATVALUE(LINDATA,CODE,IDX)
%
%CODE: field name in the LINDATA structure. May be:
%
%   ElemIndex   - ordinal position in the RING 
%   SPos        - longitudinal position [m]
%   ClosedOrbit - closed orbit column vector with 
%                 components x, px, y, py (momentums, NOT angles)						
%   Dispersion  - dispersion orbit position vector with 
%                 components eta_x, eta_prime_x, eta_y, eta_prime_y
%                 calculated with respect to the closed orbit with 
%                 momentum deviation DP
%   gamma       - gamma parameter of the transformation to eigenmodes 
%   mu          - [ mux, muy] horizontal and vertical betatron phase
%   beta        - [betax, betay] vector
%   alpha       - [alphax, alphay] vector
%   emit66	- 6x6 emittance projections on x and y + energy spread
%   emit44	- emittances of the projections of beam44 on x and y
%   modemit	- [emitA emitB] emittance of modes A and B (should be constant)

%IDX:   index of desired values (default: all)

switch code
   case {'ClosedOrbit','Dispersion'}
   v=cat(2,lindata.(code))';
   otherwise
   v=cat(1,lindata.(code));
end
if nargin > 2
   v=v(:,idx);
end
