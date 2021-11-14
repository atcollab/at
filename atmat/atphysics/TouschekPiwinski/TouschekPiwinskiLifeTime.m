function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(r,dpp,Ib,varargin)
% function [Tl,contributionsTL]=TouschekPiwinskiLifeTime(ring,dpp,Ib,...)
%
% evaluates Touschek Lifetime using Piwinski formula
%
% INPUT
%
% TouschekPiwinskiLifeTime(
%  ring,
%  momentumaperturecolumnvector,  column array (size of r or positions)
%                                 it can be length(r)x2, for positive and
%                                 negative aperture
%  current per bunch in A,        scalar
%  positions where to evaluate,	  default: all elements with length>0  column array
%  emittancex,                    default: atx modemittance(1)   scalar
%  emittancey,                    default: emittancex/2		     scalar
%  integration_method,            default: 'integral', may be: 'integral', 'quad', 'trapz')
%  energy_spread,                 scalar
%  bunch_length,	              scalar
%			   )
%
% OUTPUT
%
%  contributionsTL 1/T contribution at each element 
%                  if dpp has positive and negative apertures, then contributionTL is a 2 columns vector      
%
%  Tl  Lifetime in seconds 1/Tl=sum(contributionsTL.*L)/sum(L);
%
% NAME-VALUE PAIRS
%
% TouschekPiwinskiLifeTime(..., 'AbsTol', abstol)
%   Absolute tolerance for the 'integral' function. Default: 1.0e-16
%
% TouschekPiwinskiLifeTime(..., 'RelTol', abstol)
%   Relative tolerance for the 'integral' function. Default: 1.0e-12
%
% "The Touscheck Effect in strong focusing storage rings"
% A.Piwinski, DESY 98-179, November 1998
%
% "Touscheck Effect calculation and its applications to a transport line"
% A.Xiao M. Borland, Proceedings of PAC07, Albuquerque, New Mexico, USA
%

% created 31-10-2012
% updated 22-01-2013 corrected and compared with elegant
% updates 08-11-2018 add: PhysConstant, default 'integral', contributionTL
%                         given for positive and negative side
% updated 14/11/2021 Added optional specification of integral tolerances

%ensure a column lattice
r=reshape(r,numel(r),1);

e0 = PhysConstant.elementary_charge.value; %1.60217646e-19; %Coulomb
r0 = PhysConstant.classical_electron_radius.value; %2.817940327e-15; %m %  classical electron radius
spl = PhysConstant.speed_of_light_in_vacuum.value; %299792458; % speed of ligth

[abstol,varargs]=getoption(varargin, 'AbsTol', 1.0e-16);
[reltol,varargs]=getoption(varargs, 'RelTol', 1.0e-12);
tol = {'AbsTol', abstol, 'RelTol', reltol};

[positions,varargs]=getargs(varargs,[]);
if isempty(positions)
    % positions default= non zero length elements
    positions=findcells(r,'Length');
    L=getcellstruct(r,'Length',positions);
    positions=positions(L>0);
    size(positions);
end
% get optics
[lo,pa]=atx(r,0,positions);

emitx=pa.modemittance(1);
emity=emitx./2;

[emitx,emity,integrationmethod,sigp,sigs,rest] = ...
    getargs(varargs, emitx, emity,'integral',pa.espread,pa.blength);

fprintf('emitx: %.3e [m]\n', emitx);
fprintf('emity: %.3e [m]\n', emity);
fprintf('integration method: "%s"\n', integrationmethod);
fprintf('energy spread: %.3e\n', sigp);
fprintf('bunch length:  %.5g [m]\n', sigs);

% if naddvar==2
%     emitx=varargs{2};
%     disp('set defaults: ey=ex/2')
%     disp(' integration method is integral,')
%     disp(' energy spread, bunch length from ATX')
%     
% elseif naddvar==3
%     emitx=varargs{2};
%     emity=varargs{3};
%     disp('set defaults: ')
%     disp(' integration method is integral,')
%     disp(' energy spread, bunch length from ATX')
%     
% elseif naddvar==4
%     emitx=varargs{2};
%     emity=varargs{3};
%     integrationmethod=varargs{4};
%     disp('set defaults: ')
%     disp(' energy spread, bunch length from ATX')
%     
% elseif naddvar==5
%     emitx=varargs{2};
%     emity=varargs{3};
%     integrationmethod=varargs{4};
%     sigp=varargs{5};
%     disp('set defaults: ')
%     disp('bunch length from ATX')
%     
% elseif naddvar==6
%     emitx=varargs{2};
%     emity=varargs{3};
%     integrationmethod=varargs{4};
%     sigp=varargs{5};
%     sigs=varargs{6};
%     
% else
%     disp('set defaults: ey=ex/2')
%     disp(' integration method is integral,')
%     disp(' energy spread, bunch length, x emittance from ATX')
%     disp(' evaluation at all points with non zero length')
% end

% if dpp is a scalar assume constant momentum aperture.
if numel(dpp)==1
    dpp=dpp*ones(size(positions'));
end

dppinput=dpp;
Tlcol=zeros(1,size(dppinput,2));
Circumference=findspos(r,length(r)+1);

E0=atenergy(r);
Nb = Ib/(spl/Circumference)/e0; %Number of particle per bunch.
mass_el_ev=PhysConstant.electron_mass_energy_equivalent_in_MeV.value*1e6;
relgamma = E0/mass_el_ev;%0.510998928e6;
relbeta=sqrt(1-1./relgamma.^2);

aaa=cat(1,lo.alpha);
bbb=cat(1,lo.beta);
ddd=cat(2,lo.Dispersion)';

bx=bbb(:,1); % betax
by=bbb(:,2); % betay
Dx=ddd(:,1);
Dy=ddd(:,3);

ax=aaa(:,1);
ay=aaa(:,2);
Dpx=ddd(:,2);
Dpy=ddd(:,4);

sigxb=sqrt(emitx.*bx);
sigyb=sqrt(emity.*by);

sigx=sqrt(emitx.*bx+sigp.^2.*Dx.^2);
sigy=sqrt(emity.*by+sigp.^2.*Dy.^2); %%%  was mistaken! it was Dx!!!!!!

Dtx=Dx.*ax+Dpx.*bx;%  % alpha=-b'/2
Dty=Dy.*ay+Dpy.*by;%

%     sigtx=sqrt(sigx.^2+sigp.^2.*Dtx.^2);
%     sigty=sqrt(sigy.^2+sigp.^2.*Dty.^2);
% sigtx2=sigx.^2+sigp.^2.*Dtx.^2;
%     sigty2=sigy.^2+sigp.^2.*Dty.^2;

sigp2=sigp.^2;
Dx2=Dx.^2;
Dy2=Dy.^2;
Dtx2=Dtx.^2;
Dty2=Dty.^2;
sigxb2=sigxb.^2;
sigyb2=sigyb.^2;

sighinv2=1./(sigp2) +(Dx2+Dtx2)./(sigxb2) + (Dy2+Dty2)./(sigyb2);
sigh=sqrt(1./sighinv2);

B1=1./(2.*(relbeta.^2).*(relgamma.^2)).*( (bx.^2./(sigxb.^2)).*(1-(sigh.^2.*Dtx.^2./(sigxb.^2))) + (by.^2./(sigyb.^2)).*(1-(sigh.^2.*Dty.^2./(sigyb.^2))));
B2sq=1./(4.*(relbeta.^4).*(relgamma.^4)).*((bx.^2./(sigxb.^2)).*(1-(sigh.^2.*Dtx.^2./(sigxb.^2)))-(by.^2./(sigyb.^2)).*(1-(sigh.^2.*Dty.^2./(sigyb.^2)))).^2+(sigh.^4.*bx.^2.*by.^2.*Dtx.^2.*Dty.^2)./((relbeta.^4).*(relgamma.^4).*sigxb.^4.*sigyb.^4);
B2=sqrt(B2sq);
contributionsTL=zeros(size(dppinput));

for dppcolnum=1:size(dppinput,2)
    dpp=dppinput(:,dppcolnum);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%
    %%%%%%%% From here calculation takes place.
    %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    um=relbeta.^2*dpp.^2;
    %  em=bx.^2.*sigx.^2./(relbeta.^2.*relgamma.^2.*sigxb.^2.*sigtx2).*um;
    
    val=zeros(size(B1));
    km=atan(sqrt(um));
    FpiWfact=sqrt(pi.*(B1.^2-B2.^2)).*um;
    
    for ii=1:length(positions)
        % choose integration method
        switch integrationmethod
            case 'infiniteintegral'
                val(ii)= integral(@(u)TLT_IntPiw(u,um(ii),B1(ii),B2(ii)),um(ii),Inf, tol{:}); %,...um(ii)*1e4
                
            case 'integral'
                val(ii) = integral(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2, tol{:}); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            case 'quad'
                val(ii)= quad(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            case 'trapz'
                k=linspace(km(ii),pi/2,10000);
                val(ii)= trapz(k,TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii))); %,...,'Waypoints',evalpoints um(ii)*1e4
                
            otherwise % use default method quad (backward compatible)
                val(ii)=integral(@(k)TLT_IntPiw_k(k,km(ii),B1(ii),B2(ii)),km(ii),pi/2, tol{:}); %,...,'Waypoints',evalpoints um(ii)*1e4
                
        end
    end
    
    
    
    switch integrationmethod
        case 'infiniteintegral'
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*FpiWfact;
            
        case {'integral' 'quad' }
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*2.*FpiWfact;
        case {'trapz'}
            frontfact=(r0.^2.*spl.*Nb.*sigh.*bx.*by)./(4.*sqrt(pi).*(relbeta.^2).*(relgamma.^4).*sigxb.^2.*sigyb.^2.*sigs.*sigp);
            
        otherwise
            frontfact=(r0.^2.*spl.*Nb)./(8.*pi.*(relgamma.^2).*sigs.*sqrt(...
                (sigx.^2).*(sigy.^2)-sigp.^4.*Dx.^2.*Dy.^2).*um).*2.*FpiWfact;
            
    end
    contributionsTL(:,dppcolnum)=frontfact.*val;
    
    L=getcellstruct(r,'Length',positions);
    Tlcol(dppcolnum)=1/(1/sum(L)*sum(contributionsTL(:,dppcolnum).*L));
    
end

Tl=length(Tlcol)/(sum(1./Tlcol));

return
