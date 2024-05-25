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
    initial=[long(needed(1:end-1));false]; %needed
    final=[false;initial(1:end-1)]; %needed

    lg=initial(refs(needed)); % refpts
    L=atgetfieldvalues(ring(long),'Length'); %long
    
    beta0=avebeta(lg,:); %long
    alpha0=cat(1,lind(initial).alpha); %long
    mu0=avemu(lg,:); %long
    disp0=avedisp(lg,:); %long
    ClosedOrbit0=ClosedOrbit(lg,:); %long
    
    mu1=cat(1,lind(final).mu); %long
    disp1=cat(2,lind(final).Dispersion)'; %long
    ClosedOrbit1=cat(2,lind(final).ClosedOrbit)'; %long
    
    L2=[L L]; %long
    avebeta(lg,:)=betadrift(beta0,alpha0,L2);
    avemu(lg,:)=0.5*(mu0+mu1);
    avedisp(lg,[1 3])=(disp1(:,[1 3])+disp0(:,[1 3]))*0.5;
    avedisp(lg,[2 4])=(disp1(:,[1 3])-disp0(:,[1 3]))./L2;
    foc=atgetcells(ring(long),'PolynomB',@(el,polb) length(polb)>=2 && polb(2)~=0||polb(3)~=0); %long
    if any(foc)
        qp=false(size(lg));
        qp(lg)=foc;

        %Extract element parameters
        reng_selection=ring(refpts(qp));
        L=atgetfieldvalues(reng_selection,'Length','Default',0);
        q=eps()+atgetfieldvalues(reng_selection,'PolynomB',{2},'Default',0);
        m=atgetfieldvalues(reng_selection,'PolynomB',{3},'Default',0);
        R11=atgetfieldvalues(reng_selection,'R2',{1,1},'Default',1);
        dx=(atgetfieldvalues(reng_selection,'T2',{1},'Default',0)-atgetfieldvalues(reng_selection,'T1',{1},'Default',0))/2;
        ba=atgetfieldvalues(reng_selection,'BendingAngle','Default',0);
        irho=ba./L;
        e1=atgetfieldvalues(reng_selection,'EntranceAngle','Default',0);
        Fint=atgetfieldvalues(reng_selection,'Fint','Default',0);
        gap=atgetfieldvalues(reng_selection,'gap','Default',0);
        
        %Hard edge model on dipoles
        d_csi=ba.*gap.*Fint.*(1+sin(e1).^2)./cos(e1)./L;
        Cp=[ba.*tan(e1)./L -ba.*tan(e1-d_csi)./L];
        for ii=1:2
            alpha0(foc,ii)=alpha0(foc,ii)-beta0(foc,ii).*Cp(:,ii);
            disp0(foc,2+(ii-1)*2)=disp0(foc,2+(ii-1)*2)-disp0(foc,1+(ii-1)*2).*Cp(:,ii);
        end
        
        % Additional components
        dx0=(ClosedOrbit0(foc,1)+ClosedOrbit1(foc,1))/2;
        K=q.*R11+2*m.*(dx0-dx);
        K2=[K -K]; %long
        K2(:,1)=K2(:,1)+irho.^2;
        irho2=[irho 0*irho];

        % Apply formulas
        avebeta(qp,:)=betafoc(beta0(foc,:),alpha0(foc,:),K2,L2(foc,:));
        avedisp(qp,:)=dispfoc(disp0(foc,:),irho2,K2,L2(foc,:));
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
            disp0(:,[2 4]).*(1-cos(sqrt(K).*L))./K)./L+...
            ir.*(L-sin(sqrt(K).*L)./sqrt(K))./K./L;
        avedisp(:,[2 4])=(disp0(:,[2 3]).*(sin(sqrt(K).*L)./sqrt(K))-...
            disp0(:,[1 3]).*(1-cos(sqrt(K).*L)))./L+...
            ir.*(L-cos(sqrt(K).*L))./K./L;
    end
end
