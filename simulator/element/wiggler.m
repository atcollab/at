function [z] = wiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method)
% wiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method)
%
% FamName	family name
% Ltot		total length of the wiggle
% Lw		total length of the wiggle
% Bmax	 	peak wiggler field [Tesla]
% Nstep		num of integration steps per period
% Nmeth		symplectic integration method, 2nd or 4th order: 2 or 4
% By		wiggler harmonics for horizontal wigglers
% Bx		wiggler harmonics for vertical wigglers
% method        name of the function to use for tracking
%
% returns assigned address in the FAMLIST that is uniquely identifies
% the family

%---------------------------------------------------------------------------
% Modification Log:
% -----------------
% .03  2003-06-19	YK Wu, Duke University, wu@fel.duke.edu
%                               Add checks for input arguments
% .02  2003-06-18	YK Wu, Duke University
%				Add checks for inputs, add comments
%
% .01  2003-04-20	YK Wu, J. Li, Duke University  
%				Define a wiggler element
%
%---------------------------------------------------------------------------
%  Accelerator Physics Group, Duke FEL Lab, www.fel.duke.edu
%

global MaxOrder;
global NumIntSteps;

GWIG_EPS = 1e-6;
dNw = abs(mod(Ltot/Lw, 1));
if dNw > GWIG_EPS
  error(' Wiggler: Ltot/Lw is not an integter.');
end

ElemData.FamName        = fname;  % add check for identical family names
ElemData.Length		= Ltot;
ElemData.Lw             = Lw;
ElemData.Bmax           = Bmax;
ElemData.Nstep    	= Nstep;
ElemData.Nmeth      	= Nmeth;
if ~isempty(By)
  ElemData.NHharm       = length(By(1,:));
  for i=1:ElemData.NHharm
    kx = By(3,i); ky = By(4,i); kz = By(5,i);
    dk = sqrt(abs(ky*ky - kz*kz - kx*kx))/abs(kz);
    if ( dk > GWIG_EPS ) then
      error([' Wiggler (H): kx^2 + kz^2 - ky^2 != 0!, i = ', num2str(i,3)]);
    end;
  end
else
  ElemData.NHharm         = 0;
end

if ~isempty(Bx)
  ElemData.NVharm         = length(Bx(1,:));
  for i=1:ElemData.NVharm
    kx = Bx(3,i); ky = Bx(4,i); kz = Bx(5,i);
    dk = sqrt(abs(kx*kx - kz*kz - ky*ky))/abs(kz);
    if ( dk > GWIG_EPS ) then
      error([' Wiggler (V): ky^2 + kz^2 - kx^2 != 0!, i = ', num2str(i,3)]);
    end;
  end
else
  ElemData.NVharm         = 0;
end
ElemData.By             = By;
ElemData.Bx             = Bx;
ElemData.R1             = diag(ones(6,1));
ElemData.R2             = diag(ones(6,1));
ElemData.T1             = zeros(1,6);
ElemData.T2             = zeros(1,6);
ElemData.PassMethod 	= method;


global FAMLIST
z = length(FAMLIST)+1; % number of declare families including this one
FAMLIST{z}.FamName = fname;
FAMLIST{z}.NumKids = 0;
FAMLIST{z}.KidsList= [];
FAMLIST{z}.ElemData= ElemData;
