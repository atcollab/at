function [lindata,avebeta,avemu,avedisp,nu,xsi]=atavedata(ring,dpp,refpts,varargin)
%ATAVEDATA       Average of optical functions on selected elements
%
%[LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS)
%
%LINDATA : Identical to ATLINOPT output
%AVEBEA :  Average Beta functions
%AVEMU :   Average phase advance
%AVEDISP : Average dispersion
%TUNES : Vector of tunes
%CHROMS : Vector of chromaticites
%
%[LINDATA,AVEBETA,AVEMU,AVEDISP,TUNES,CHROMS]=ATAVEDATA(RING,DPP,REFPTS,ORBITIN)
%    does not search for closed orbit. instead ORBITIN is used


lr=length(ring)+1;
if islogical(refpts)
    refs=[refpts(:);false(lr-length(refpts),1)];
else
    refs=false(lr,1); % lr
    refs(refpts)=true;
end

long=atgetcells(ring,'Length',@(elem,lg) lg>0) & refs(1:end-1); %lr-1
needed=refs | [false;long]; %lr
initial=[long(needed(1:end-1));false]; %needed
final=[false;initial(1:end-1)]; %needed
setoption('WarningDp6D',false)
[ringdata, lind]=atlinopt6(ring,needed,'dp',dpp,'get_chrom',varargin{:}); %needed
setoption('WarningDp6D',true)
nu = ringdata.tune;
xsi = ringdata.chromaticity;

lindata=lind(refs(needed)); %refpts
avebeta=cat(1,lindata.beta); %refpts
avemu=cat(1,lindata.mu); %refpts
avedisp=cat(2,lindata.Dispersion)'; %refpts

ringlong = ring(long);
lin0 = lind(initial);
lin1 = lind(final);

lg=initial(refs(needed)); % refpts

beta0=avebeta(lg,:); %long
alpha0=cat(1,lin0.alpha); %long
mu0=avemu(lg,:); %long
disp0=avedisp(lg,:); %long

mu1=cat(1,lin1.mu); %long
ClosedOrbit0=cat(2,lin0.ClosedOrbit)'; %long
ClosedOrbit1=cat(2,lin1.ClosedOrbit)'; %long
dx0=(ClosedOrbit0(:,1)+ClosedOrbit1(:,1))/2;

[L,q,m,ba,R11,dx] = getparams(ringlong, @longparams);
irho=ba./L;
K = q.*R11 + 2*m.*(dx0-dx);
Kx = K+irho.*irho;
Ky = -K;
bend = (irho ~= 0.0);

% Hard edge model on dipoles
[e1,Fint,gap] = getparams(ringlong(bend), @focparams);
d_csi=irho(bend).*gap.*Fint.*(1+sin(e1).^2)./cos(e1);
Cp=[irho(bend).*tan(e1) -irho(bend).*tan(e1-d_csi)];
alpha0(bend,:)=alpha0(bend,:)-beta0(bend,:).*Cp;
disp0(bend,[2 4])=disp0(bend,[2 4])-disp0(bend,[1 3]).*Cp;

avemu(lg,:)=0.5*(mu0+mu1);
avebeta(lg,:) = betalong(beta0,alpha0,[Kx Ky],[L L]);
avedisp(lg,1:2)=displong(disp0(:,1:2),irho,Kx,L);
avedisp(lg,3:4)=displong(disp0(:,3:4),zeros(size(irho)),Ky,L);

    function avebeta=betadrift(beta0,alpha0,L)
        gamma0=(alpha0.*alpha0+1)./beta0;
        avebeta=beta0-alpha0.*L+gamma0.*L.*L/3;
    end

    function avebeta=betafoc(beta0,alpha0,K,L)
        gamma0=(alpha0.*alpha0+1)./beta0;
        avebeta=((beta0+gamma0./K).*L+(beta0-gamma0./K).*sin(2.0*sqrt(K).*L)./sqrt(K)/2.0+...
            (cos(2.0*sqrt(K).*L)-1.0).*alpha0./K)./L/2.0;
    end

    function avebeta=betalong(beta0,alpha0,K,L)
        avebeta=zeros(size(beta0));
        kp = abs(K) >= 1.0e-7;
        k0 = ~kp;
        avebeta(kp) = betafoc(beta0(kp),alpha0(kp),K(kp),L(kp));
        avebeta(k0) = betadrift(beta0(k0),alpha0(k0),L(k0));
    end


    function avedisp=dispfoc(disp0,ir,K,L)
        eta=(disp0(:,1).*(sin(sqrt(K).*L)./sqrt(K))+...
            disp0(:,2).*(1-cos(sqrt(K).*L))./K+...
            ir.*(L-sin(sqrt(K).*L)./sqrt(K))./K)./L;
        etap=(disp0(:,2).*(sin(sqrt(K).*L)./sqrt(K))-...
            disp0(:,1).*(1-cos(sqrt(K).*L))+...
            ir.*(1-cos(sqrt(K).*L))./K)./L;
        avedisp=[eta etap];
    end

    function avedisp=dispdrift(disp0,ir,L)
        eta=disp0(:,1)+disp0(:,2).*L/2+ir.*L.*L./6.0;
        etap=disp0(:,2)+ir.*L/2.0;
        avedisp=[eta etap];
    end

    function avedisp=displong(disp0,ir,K,L)
        avedisp=zeros(size(disp0));
        kp = abs(K) >= 1.0e-7;
        k0 = ~kp;
        avedisp(kp,:) = dispfoc(disp0(kp,:),ir(kp),K(kp),L(kp));
        avedisp(k0,:) = dispdrift(disp0(k0,:),ir(k0),L(k0));
    end

    function [L,q,m,ba,R11,dx] = longparams(elem)
        m = 0;
        q = 0;
        if isfield(elem, 'PolynomB')
            pb = elem.PolynomB;
            if length(pb) >= 3
                m = pb(3);
                q = pb(2);
            elseif length(pb) >= 2
                q = pb(2);
            end
        end
        if isfield(elem, 'Length'), L = elem.Length; else, L = 0; end
        if isfield(elem, 'BendingAngle'), ba = elem.BendingAngle; else, ba = 0; end
        if isfield(elem, 'R2'), R11 = elem.R2(1,1); else, R11 = 1; end
        if isfield(elem, 'T1'), x1 = elem.T1(1); else, x1 = 0; end
        if isfield(elem, 'T2'), x2 = elem.T2(1); else, x2 = 0; end
        dx = (x2 - x1)/2;
    end

    function [e1,Fint, gap] = focparams(elem)
        if isfield(elem, 'EntranceAngle'), e1 = elem.EntranceAngle; else, e1 = 0; end
        if isfield(elem, 'Fint'), Fint = elem.Fint; else, Fint = 0; end
        if isfield(elem, 'gap'), gap = elem.gap; else, gap = 0; end
    end

    function varargout = getparams(ring, func)
        [varargout{1:nargout}] = cellfun(func, ring);
    end

end
