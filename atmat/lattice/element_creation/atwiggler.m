function elem=atwiggler(fname,varargin)
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
% .04  2018-07-30   A.Mash'al, Iranian Light Source Facility 
%                               Add energy to ElemData
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
[rsrc,Ltot,Lw,Bmax,Nstep,Nmeth,By,Bx,method] = decodeatargs({0,0,0,0,0,zeros(6,1),zeros(6,1),'WiggLinearPass'},varargin);
[Ltot,rsrc] = getoption(rsrc,'Ltot',Ltot);
[Lw,rsrc] = getoption(rsrc,'Lw',Lw);
[Bmax,rsrc] = getoption(rsrc,'Bmax',Bmax);
[Nstep,rsrc] = getoption(rsrc,'Nstep',Nstep);
[Nmeth,rsrc] = getoption(rsrc,'Nmeth',Nmeth);
[Bx,rsrc] = getoption(rsrc,'Bx',Bx);
[By,rsrc] = getoption(rsrc,'By',By);
[method,rsrc] = getoption(rsrc,'PassMethod',method);
[cl,rsrc] = getoption(rsrc,'Class','Wiggler');
GWIG_EPS = 1e-6;
dNw = abs(mod(Ltot/Lw, 1));
if dNw > GWIG_EPS
  error(' Wiggler: Ltot/Lw is not an integer.');
end

if ~isempty(By)
  NHharm = length(By(1,:));
  for i=1:NHharm
    kx = By(3,i); ky = By(4,i); kz = By(5,i);
    dk = sqrt(abs(ky*ky - kz*kz - kx*kx))/abs(kz);
    if ( dk > GWIG_EPS ) 
      error([' Wiggler (H): kx^2 + kz^2 - ky^2 != 0!, i = ', num2str(i,3)]);
    end
  end
else
  NHharm = 0;
end

if ~isempty(Bx)
  NVharm = length(Bx(1,:));
  for i=1:NVharm
    kx = Bx(3,i); ky = Bx(4,i); kz = Bx(5,i);
    dk = sqrt(abs(kx*kx - kz*kz - ky*ky))/abs(kz);
    if ( dk > GWIG_EPS ) 
      error([' Wiggler (V): ky^2 + kz^2 - kx^2 != 0!, i = ', num2str(i,3)]);
    end
  end
else
  NVharm = 0;
end


elem = atbaselem(fname,method,'Class',cl,'Ltot',Ltot,'Lw',Lw,...
    'Bmax',Bmax,'Nstep',Nstep,'Nmeth',Nmeth,'Bx',Bx,'By',By,...
    'NHharm',NHharm','NVharm',NVharm,rsrc{:});
end
