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
[lind,nu,xsi]=atlinopt(ring,dpp,find(needed),varargin{:}); %needed

lindata=lind(refs(needed)); %refpts
avebeta=cat(1,lindata.beta); %refpts
avemu=cat(1,lindata.mu); %refpts
avedisp=cat(2,lindata.Dispersion)'; %refpts
ClosedOrbit=cat(2,lindata.ClosedOrbit)'; %refpts


if any(long)
    ringlong = ring(long);
    initial=[long(needed(1:end-1));false]; %needed
    final=[false;initial(1:end-1)]; %needed

    lg=initial(refs(needed)); % refpts
    L=atgetfieldvalues(ringlong,'Length'); %long
    
    beta0=avebeta(lg,:); %long
    alpha0=cat(1,lind(initial).alpha); %long
    mu0=avemu(lg,:); %long
    disp0=avedisp(lg,:); %long
    ClosedOrbit0=ClosedOrbit(lg,:); %long
    
    mu1=cat(1,lind(final).mu); %long
    disp1=cat(2,lind(final).Dispersion)'; %long
    ClosedOrbit1=cat(2,lind(final).ClosedOrbit)'; %long
    dx0=(ClosedOrbit0(:,1)+ClosedOrbit1(:,1))/2;
    
    L2=[L L]; %long
    avebeta(lg,:)=betadrift(beta0,alpha0,L2);
    avemu(lg,:)=0.5*(mu0+mu1);
    avedisp(lg,[1 3])=(disp1(:,[1 3])+disp0(:,[1 3]))*0.5;
    avedisp(lg,[2 4])=(disp1(:,[1 3])-disp0(:,[1 3]))./L2;
    L=atgetfieldvalues(ringlong,'Length','Default',0);
    q=atgetfieldvalues(ringlong,'PolynomB',{2},'Default',0);
    m=atgetfieldvalues(ringlong,'PolynomB',{3},'Default',0);
    R11=atgetfieldvalues(ringlong,'R2',{1,1},'Default',1);
    dx=(atgetfieldvalues(ringlong,'T2',{1},'Default',0)-atgetfieldvalues(ringlong,'T1',{1},'Default',0))/2;
    ba=atgetfieldvalues(ringlong,'BendingAngle','Default',0);
    irho=ba./L;
    e1=atgetfieldvalues(ringlong,'EntranceAngle','Default',0);
    Fint=atgetfieldvalues(ringlong,'Fint','Default',0);
    gap=atgetfieldvalues(ringlong,'gap','Default',0);
    K = q.*R11 + 2*m.*(dx0-dx);
    foc = abs(K) > 1.e-7;
    %Hard edge model on dipoles
    % d_csi=ba.*gap.*Fint.*(1+sin(e1).^2)./cos(e1)./L;
    % Cp=[ba.*tan(e1)./L -ba.*tan(e1-d_csi)./L];
    % alpha0=alpha0-beta0.*Cp;
    % for ii=1:2
    %     disp0(:,2*ii)=disp0(:,2*ii)-disp0(:,2*ii-1).*Cp(:,ii);
    % end

    if any(foc)
        rqp=lg;
        rqp(lg)=foc;
        kfoc = K(foc);
        irhofoc = irho(foc);
        K2 = [kfoc+irhofoc.*irhofoc -kfoc];
        irho2 = [irhofoc zeros(size(irhofoc))];
        % Apply formulas
        avebeta(rqp,:)=betafoc(beta0(foc,:),alpha0(foc,:),K2,L2(foc,:));
        avedisp(rqp,:)=dispfoc(disp0(foc,:),irho2,K2,L2(foc,:));
    end
    
end

    function avebeta=betadrift(beta0,alpha0,L)
        gamma0=(alpha0.*alpha0+1)./beta0;
        avebeta=beta0-alpha0.*L+gamma0.*L.*L/3;
    end

    function avebeta=betafoc(beta0,alpha0,K,L)
        gamma0=(alpha0.*alpha0+1)./beta0;
        avebeta=((beta0+gamma0./K).*L+(beta0-gamma0./K).*sin(2.0*sqrt(K).*L)./sqrt(K)/2.0+...
            (cos(2.0*sqrt(K).*L)-1.0).*alpha0./K)./L/2.0;
    end

    function avedisp=dispfoc(disp0,ir,K,L)
        avedisp=disp0;
        avedisp(:,[1 3])=(disp0(:,[1 3]).*(sin(sqrt(K).*L)./sqrt(K))+...
            disp0(:,[2 4]).*(1-cos(sqrt(K).*L))./K+...
            ir.*(L-sin(sqrt(K).*L)./sqrt(K))./K)./L;
        avedisp(:,[2 4])=(disp0(:,[2 4]).*(sin(sqrt(K).*L)./sqrt(K))-...
            disp0(:,[1 3]).*(1-cos(sqrt(K).*L))+...
            ir.*(1-cos(sqrt(K).*L))./K)./L;
    end
end
