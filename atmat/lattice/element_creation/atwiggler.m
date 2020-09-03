function elem=atwiggler(fname,varargin)
% ATWIGGLER Creates a wiggler
%
%ELEM=ATWIGGLER(FAMNAME, LENGTH, LW, BMAX, ENERGY, PASSMETHOD)
%
% FAMNAME       family name
% LENGTH        total length
% LW            Period length
% BMAX          Peak magnetic field [T]
% ENERGY        Beam energy [eV]
% PASSMETHOD    Tracking function. Default 'GWigSymplecticPass'
%
%ELEM=ATWIGGLER(...,'keyword',value...)
%
% Keywords:
% Nstep		number of integration steps per period (default 5)
% Nmeth		symplectic integration method, 2nd or 4th order: 2 or 4 (default 4)
% By		harmonics of the horizontal wiggler. Default [1;1;0;1;1;0]
%               6xNH matrix, with NH number of harmonics
% Bx		harmonics of the vertical wigglers. Default []
%               6xNV matrix, with NV number of harmonics
%
%see also: GWigSymplecticPass

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
[rsrc,Ltot,Lw,Bmax,energy,method] = decodeatargs({0,0,0,0,'GWigSymplecticPass'},varargin);
[Ltot,rsrc] = getoption(rsrc,'Ltot',Ltot);
[Lw,rsrc] = getoption(rsrc,'Lw',Lw);
[Bmax,rsrc] = getoption(rsrc,'Bmax',Bmax);
[Nstep,rsrc] = getoption(rsrc,'Nstep',5);
[Nmeth,rsrc] = getoption(rsrc,'Nmeth',4);
[Bx,rsrc] = getoption(rsrc,'Bx',[]);
[By,rsrc] = getoption(rsrc,'By',[1;1;0;1;1;0]);
[method,rsrc] = getoption(rsrc,'PassMethod',method);
[cl,rsrc] = getoption(rsrc,'Class','Wiggler');
GWIG_EPS = 1e-6;
dNw = Ltot/Lw-round(Ltot/Lw);
if dNw > GWIG_EPS
    error(' Wiggler: Ltot/Lw is not an integer.');
end

By=reshape(By,6,[]);
NHharm = size(By,2);
for i=1:NHharm
    kx = By(3,i); ky = By(4,i); kz = By(5,i);
    dk = sqrt(abs(ky*ky - kz*kz - kx*kx))/abs(kz);
    if ( dk > GWIG_EPS )
        error([' Wiggler (H): kx^2 + kz^2 - ky^2 != 0!, i = ', num2str(i,3)]);
    end
end

Bx=reshape(Bx,6,[]);
NVharm = size(Bx,2);
for i=1:NVharm
    kx = Bx(3,i); ky = Bx(4,i); kz = Bx(5,i);
    dk = sqrt(abs(kx*kx - kz*kz - ky*ky))/abs(kz);
    if ( dk > GWIG_EPS )
        error([' Wiggler (V): ky^2 + kz^2 - kx^2 != 0!, i = ', num2str(i,3)]);
    end
end

elem = atbaselem(fname,method,'Class',cl,'Length',Ltot,'Lw',Lw,...
    'Bmax',Bmax,'Energy',energy,'Nstep',Nstep,'Nmeth',Nmeth,'Bx',Bx,'By',By,...
    'NHharm',NHharm','NVharm',NVharm,rsrc{:});
end
