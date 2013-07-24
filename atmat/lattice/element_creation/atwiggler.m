function Elem = atwiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method)
% atwiggler(fname, Ltot, Lw, Bmax, Nstep, Nmeth, By, Bx, method)
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
% returns a wiggler structure with class 'Wiggler'

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

Elem.FamName        = fname;  % add check for identical family names
Elem.Length		= Ltot;
Elem.Lw             = Lw;
Elem.Bmax           = Bmax;
Elem.Nstep    	= Nstep;
Elem.Nmeth      	= Nmeth;
if ~isempty(By)
  Elem.NHharm       = length(By(1,:));
  for i=1:Elem.NHharm
    kx = By(3,i); ky = By(4,i); kz = By(5,i);
    dk = sqrt(abs(ky*ky - kz*kz - kx*kx))/abs(kz);
    if ( dk > GWIG_EPS ) then
      error([' Wiggler (H): kx^2 + kz^2 - ky^2 != 0!, i = ', num2str(i,3)]);
    end;
  end
else
  Elem.NHharm         = 0;
end

if ~isempty(Bx)
  Elem.NVharm         = length(Bx(1,:));
  for i=1:Elem.NVharm
    kx = Bx(3,i); ky = Bx(4,i); kz = Bx(5,i);
    dk = sqrt(abs(kx*kx - kz*kz - ky*ky))/abs(kz);
    if ( dk > GWIG_EPS ) then
      error([' Wiggler (V): ky^2 + kz^2 - kx^2 != 0!, i = ', num2str(i,3)]);
    end;
  end
else
  Elem.NVharm         = 0;
end
Elem.By             = By;
Elem.Bx             = Bx;
Elem.R1             = diag(ones(6,1));
Elem.R2             = diag(ones(6,1));
Elem.T1             = zeros(1,6);
Elem.T2             = zeros(1,6);
Elem.PassMethod 	= method;
Elem.Class          = 'Wiggler';