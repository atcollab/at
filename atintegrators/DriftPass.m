function varargout=DriftPass(varargin)
%  DriftPass - Integrator for Drift spaces
% 
%  The transformation performed is:
% 
%  x1  = L/(1+d0) * x0' + x0 
%  x1' = x0'
%  y1  = L/(1+d0) * y0' + y0 
%  y1' = y0'
%  d1  = d0
%  ct1 = L/(1+d0)^2 * (x0'^2 + y0'^2)/2 + ct0
%  
%  The list of actions performed is:
%  - translate coordinate system
%  - rotate coordinate system
%  - test aperture (rectangular, then elliptical)
%  - apply drift transformation
%  - test aperture (rectangular, then elliptical)
%  - rotate coordinate system
%  - translate coordinate system
%  
%  The calling syntax is:
% 
% 		OUTPUT_COORD = DriftPass(ELEMENT, INPUT_COORD)
%
% 		[req, opt] = DriftPass 
%   
%   where:
%   ELEMENT is a matlab structure with fields
% 
%       Required:
%       ELEMENT.Length, (double) 
%   
%       Optional:
%       ELEMENT.R1 (6x6 double)  
%       ELEMENT.R2 (6x6 double)   
%       ELEMENT.T1 (6x1 double)
%       ELEMENT.T2 (6x1 double)
%       ELEMENT.EApertures (2x1 double)
%       ELEMENT.RApertures (4x1 double)
%   
%   INPUT_COORD is a 6xN matrix of coordinates for N particles at the
%   entrance of ELEMENT
%
%   OUTPUT_COORD is a 6xN vecotr of coordinated at the exit of ELEMENT
%
%   Description of the ELEMENT fields for DriftPass
%   R1: is the coordinate rotation matrix at the entrance of ELEMENT 
%   R2: is the coordinate rotation matrix at the exit of ELEMENT 
%   T1: is the coordinate translation vector at the entrance of ELEMENT 
%   T2: is the coordinate translation vector at the exit of ELEMENT 
%   EApertures: [ax,ay] are horizontal and vertical semi axis of an
%       aperture ellipse, tested at the entrance and exit of ELEMENT
%   RApertures: [-x x -y y] are negative horizontal, positive horizontal, 
%       negative vertical and positive vertical half size of a 
%       rectangular aperture, tested at the entrance and exit of ELEMENT
%   
%   Example of use
%   >> L=0.1;
%   >> ELEMENT=atdrift('DRIFT',L);
%   >> outcoord=DriftPass(ELEMENT,[0 0 0 0 0 0]')
%
%see also: atdrift, ringpass, linepass, DriftPass.c

error('AT:MissingMexFile','"%s.%s" is missing. Did you run "atmexall" ?',...
    mfilename,mexext);
end
