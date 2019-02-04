function [...
    dX_dq, dY_dq, ...
    dDx_dq, dDx_db, Dx, ...
    dXY_ds, dYX_ds, ...
    dDx_ds, dDy_ds, ...
    dDy_da...
    ]=AnalyticResponseMatrixDerivative(...
    RING,dpp,...
    bpmlist,... % bpm
    cmlist,...  % correctors
    magnetlist,... % quadrupoles
    bendlist,...
    skewlist)
% Analytic Response Matrix Derivative for dipoles and quadrupoles
% 
% AnalyticResponseMatrixDerivative returns the derivative of the steerers
% response matrix when changing dipole angles and quadrupole gradients.
% 
% from loco_fast_ESRFtest.m by Zeus Marti (ALBA CELLS, Barcellona)
%
% https://arxiv.org/pdf/1711.06589.pdf
%

% define Params structure
qind=find(atgetcells(RING,'Class','Quadrupole'))';
nq=length(qind);
bind=bendlist;
nb=length(bind);
%bind=find(atgetcells(RING,'Class','Bend'))';
%bind=sort(bind([1:4:end,4:4:end]));

for iq=1:nq
    Params{iq}.FieldName='PolynomB';
    Params{iq}.FieldIndex={[1,2]};
    Params{iq}.ElemIndex=qind(iq);
end

for iq=1:nq
    Params{nq+iq}.FieldName='PolynomA';
    Params{nq+iq}.FieldIndex={[1,2]};
    Params{nq+iq}.ElemIndex=qind(iq)';
end

for ib=1:length(bind)
    Params{2*nq+ib}.FieldName='PolynomB';
    Params{2*nq+ib}.FieldIndex={[1,1]};
    Params{2*nq+ib}.ElemIndex=bind(ib);
end

for ib=1:length(bind)
    Params{2*nq+nb+ib}.FieldName='PolynomA';
    Params{2*nq+nb+ib}.FieldIndex={[1,1]};
    Params{2*nq+nb+ib}.ElemIndex=bind(ib);
end


nP=numel(Params);


isdipole=false(1,nP);
fitskewquad=false(1,nP);
fitskewdipole=false(1,nP);
fitquad=false(1,nP);
fitdipole=false(1,nP);
FirstIndex=zeros(1,nP);
LastIndex=zeros(1,nP);
numelements=zeros(1,nP);
K=eps+zeros(1,nP);
fields=cell(nP,1);
fieldindexs=cell(nP,1);

for ii = 1:length(Params)
    fields{ii}={Params{ii}.FieldName};
    fieldindexs{ii}={Params{ii}.FieldIndex};
    Index=[Params{ii}.ElemIndex];
    numelements(ii)=numel(unique(Index));
    FirstIndex(ii)=Index(1);
    LastIndex(ii)=Index(end);
    
    if isfield(RING{Index(1)},'BendingAngle')
        isdipole(ii)=true;
    end
    if isfield(RING{Index(1)},'PolynomB')
        K(ii)=eps+RING{Index(1)}.PolynomB(2);
    end
    if any(strcmp(fields{ii},'K'))  
        fitquad(ii)=true;
    end
    isB=strcmp(fields{ii},'PolynomB');
    isA=strcmp(fields{ii},'PolynomA');
    
    if any(strcmp(fields{ii},'sK'))  
        fitskewquad(ii)=true;
    elseif any(isB)
        inds=cell2mat([fieldindexs{ii}{isB}]);
        if inds(2)==1
            fitdipole(ii)=true;
        elseif inds(2)==2
            fitquad(ii)=true;
        else
            error(' fieldindexs above 2 are not supported yet... ')
        end
    elseif any(isA)
        inds=cell2mat([fieldindexs{ii}{isA}]);
        if inds(2)==1
            fitskewdipole(ii)=true;
        elseif inds(2)==2
            fitskewquad(ii)=true;
        else
            error(' fieldindexs above 2 are not supported yet... ')
        end
    else
        error('\nERROR! %s field not implemented yet! Fit Parameter %d not asigned!\n???',fields{ii}{1},ii);
    end
end

L=numelements.*getcellstruct(RING,'Length',FirstIndex)'; %%%!!! WARNING ASSUMES LENGTH EQUAL FOR ALL IN LIST?
ang=0*L;
%N.Params
Ld=2*getcellstruct(RING,'Length',bendlist)';
Kd=getcellstruct(RING,'PolynomB',bendlist,1,2)';
gap=getcellstruct(RING,'FullGap',bendlist,1,1)';
e1=getcellstruct(RING,'EntranceAngle',bendlist,1,1)';
Fint=getcellstruct(RING,'FringeInt1',bendlist,1,1)';
angd=2*getcellstruct(RING,'BendingAngle',bendlist,1,1)';
    
if isnan(Fint)
    
%    warning('no FringeInt1 defined');
    Fint=zeros(size(Fint));

end

% optics at magnet entrance (better than averge optics below)
[l,t,c]=atlinopt(RING,dpp,1:length(RING)+1);
gt.phix=arrayfun(@(a)a.mu(1),l)'/2/pi;
gt.phiy=arrayfun(@(a)a.mu(2),l)'/2/pi;
gt.betax=arrayfun(@(a)a.beta(1),l)';
gt.alfax=arrayfun(@(a)a.alpha(1),l)';
gt.etax=arrayfun(@(a)a.Dispersion(1),l)';
gt.etapx=arrayfun(@(a)a.Dispersion(2),l)';
gt.gammax=(1+gt.alfax.^2)./gt.betax; % wiedeman eq 5.98
gt.betay=arrayfun(@(a)a.beta(2),l)';
gt.alfay=arrayfun(@(a)a.alpha(2),l)';
gt.gammay=(1+gt.alfay.^2)./gt.betay;
gt.etay=arrayfun(@(a)a.Dispersion(3),l)';
gt.etapy=arrayfun(@(a)a.Dispersion(4),l)';

% % average optics inside elements
% [l,t,c]=atlinopt(RING,0,1:length(RING)+1);
% [~,avebeta,avemu,avedisp,~,~]=atavedata_mod(RING,0,(1:length(RING)+1)');
% gt.phix=avemu(:,1)/2/pi;
% gt.phiy=avemu(:,2)/2/pi;
% gt.betax=avebeta(:,1);
% gt.alfax=arrayfun(@(a)a.alpha(1),l)';
% gt.etax=avedisp(:,1);
% gt.etapx=arrayfun(@(a)a.Dispersion(2),l)';
% gt.gammax=(1+gta.alfax.^2)./gta.betax; % wiedeman eq 5.98
% gt.betay=avebeta(:,2);
% gt.alfay=arrayfun(@(a)a.alpha(2),l)';
% gt.gammay=(1+gta.alfay.^2)./gta.betay;
% gt.etay=arrayfun(@(a)a.Dispersion(3),l)';
% gt.etapy=arrayfun(@(a)a.Dispersion(4),l)';

% optics at lattice elements center
machx2=cellfun(@(a)atdivelem(a,[0.5,0.5])',RING,'un',0); machx2=[machx2{:}, atmarker('end')]';
[lx2,tx2,cx2]=atlinopt(machx2,dpp,1:length(machx2)+1);
gtc.phix=arrayfun(@(a)a.mu(1),lx2(2:2:length(machx2)+1))'/2/pi;
gtc.betax=arrayfun(@(a)a.beta(1),lx2(2:2:length(machx2)+1))';

% % average optics over elements
% [~,avebeta,avemu,avedisp,~,~]=atavedata_mod(RING,dpp,(1:length(RING)+1)');
% gtc.phix=avemu(:,1)/2/pi;
% gtc.betax=avebeta(:,1);


Qx=gt.phix(end);
Qy=gt.phiy(end);


sx=sin(pi*Qx);
cx=cos(pi*Qx);
sx2=sin(2*pi*Qx);
sy=sin(pi*Qy);
cy=cos(pi*Qy);
sy2=sin(2*pi*Qy);

% added for debug!
isdipole=fitdipole;

isdipolelist=FirstIndex(isdipole);
fitdipolelist=FirstIndex(fitdipole);
fitquadlist=FirstIndex(fitquad);
fitskewdipolelist=FirstIndex(fitskewdipole);
fitskewquadlist=FirstIndex(fitskewquad);

%bendlist=isdipolelist;

nb=numel(bpmlist);
nc=numel(cmlist);
nbd=numel(bendlist);

bxd=gt.betax(bendlist)';
axd0=gt.alfax(bendlist)';
gxd0=gt.gammax(bendlist)';
fxd=2*pi*(gt.phix(bendlist))';
dxd=gt.etax(bendlist);
dpxd=gt.etapx(bendlist);
byd=gt.betay(bendlist)';
ayd0=gt.alfay(bendlist)';
gyd0=gt.gammay(bendlist)';
fyd=2*pi*(gt.phiy(bendlist))';
dyd=gt.etay(bendlist);
dpyd=gt.etapy(bendlist);


ndq=numel(isdipolelist);
nfd=numel(fitdipolelist);
nfq=numel(fitquadlist);
nfsd=numel(fitskewdipolelist);
nfsq=numel(fitskewquadlist);

bx=gt.betax(FirstIndex)';
bxcenter=gtc.betax(FirstIndex)';
ax=gt.alfax(FirstIndex)';
gx=gt.gammax(FirstIndex)';
fx=2*pi*(gt.phix(FirstIndex))';
fxcenter=2*pi*(gtc.phix(FirstIndex))';
dx=gt.etax(FirstIndex)';
dpx=gt.etapx(FirstIndex)';
by=gt.betay(FirstIndex)';
ay=gt.alfay(FirstIndex)';
gy=gt.gammay(FirstIndex)';
fy=2*pi*(gt.phiy(FirstIndex))';
dy=gt.etay(FirstIndex)';
dpy=gt.etapy(FirstIndex)';

bxb=gt.betax(bpmlist);
fxb=2*pi*gt.phix(bpmlist);
fxbcenter=2*pi*gtc.phix(bpmlist);
byb=gt.betay(bpmlist);
fyb=2*pi*gt.phiy(bpmlist);

bxc=gt.betax(cmlist)';
fxc=2*pi*gt.phix(cmlist)';
byc=gt.betay(cmlist)';
fyc=2*pi*gt.phiy(cmlist)';

%% Some normalized parameters
Kx=K+ang.*ang./L./L;
Ky=-K;
Kdx=Kd+angd.*angd./Ld./Ld;
Kdy=-Kd;

sKx=sqrt(Kx);
sKy=sqrt(Ky);
sKlx=sqrt(Kx).*L;
sKly=sqrt(Ky).*L;

sKdx=sqrt(Kdx);
sKdy=sqrt(Kdy);
sKdlx=sqrt(Kdx).*Ld;
sKdly=sqrt(Kdy).*Ld;

%% Fringe field and edge focussing effect

d_csi=angd.*gap.*Fint.*(1+sin(e1).^2)./cos(e1)./Ld;
Cpx=angd.*tan(e1)./Ld;
Cpy=-angd.*tan(e1-d_csi)./Ld;
gxd=gxd0+bxd.*Cpx.^2-2*Cpx.*axd0;
axd=axd0-bxd.*Cpx;
gyd=gyd0+byd.*Cpy.^2-2*Cpy.*ayd0;
ayd=ayd0-byd.*Cpy;

if ndq>0
    whichbend=zeros(1,ndq);
    for jj=1:ndq
        whichbend(jj)=find(isdipolelist(jj)==bendlist);
    end
    gx(isdipole)=gxd(whichbend);
    ax(isdipole)=axd(whichbend);
    gy(isdipole)=gyd(whichbend);
    ay(isdipole)=ayd(whichbend);
    ang(isdipole)=angd(whichbend);
else
    whichbend=[];
end


%% take out quad indices
Lq=L(fitquad);
Kxq=Kx(fitquad);
Kyq=Ky(fitquad);
sKxq=sKx(fitquad);
sKyq=sKy(fitquad);
sKlxq=sKlx(fitquad);
sKlyq=sKly(fitquad);

bxq=bx(fitquad);
axq=ax(fitquad);
gxq=gx(fitquad);
fxq=fx(fitquad);
byq=by(fitquad);
ayq=ay(fitquad);
gyq=gy(fitquad);
fyq=fy(fitquad);

%% take out skew quad indices
Ls=L(fitskewquad);
if any(fitskewquad)
    Ks=zeros(1,nfsq);
    for ii=1:nfsq
        Ks(ii)=eps+RING{fitskewquadlist(ii)}.PolynomA(2);
    end
else
    Ks=Ls;
end

Kxs=Kx(fitskewquad);
Kys=Ky(fitskewquad);
sKxs=sKx(fitskewquad);
sKys=sKy(fitskewquad);
sKlxs=sKlx(fitskewquad);
sKlys=sKly(fitskewquad);

bxs=bx(fitskewquad);
axs=ax(fitskewquad);
gxs=gx(fitskewquad);
fxs=fx(fitskewquad);
dxs=dx(fitskewquad);
dpxs=dpx(fitskewquad);
bys=by(fitskewquad);
ays=ay(fitskewquad);
gys=gy(fitskewquad);
fys=fy(fitskewquad);
dys=dy(fitskewquad);
dpys=dpy(fitskewquad);


%% take out dipole indices
Lfb=L(fitdipole);
Kxfb=Kx(fitdipole);
Kyfb=Ky(fitdipole);
sKxfb=sKx(fitdipole);
sKyfb=sKy(fitdipole);
sKlxfb=sKlx(fitdipole);
sKlyfb=sKly(fitdipole);

bxfb=bx(fitdipole);
bxfbcenter=bxcenter(fitdipole);
axfb=ax(fitdipole);
gxfb=gx(fitdipole);
fxfb=fx(fitdipole);
fxfbcenter=fxcenter(fitdipole);
byfb=by(fitdipole);
ayfb=ay(fitdipole);
gyfb=gy(fitdipole);
fyfb=fy(fitdipole);

%% take out skew dipole indices
La=L(fitskewdipole);
Kxa=Kx(fitskewdipole);
Kya=Ky(fitskewdipole);
sKxa=sKx(fitskewdipole);
sKya=sKy(fitskewdipole);
sKlxa=sKlx(fitskewdipole);
sKlya=sKly(fitskewdipole);

bxa=bx(fitskewdipole);
axa=ax(fitskewdipole);
gxa=gx(fitskewdipole);
fxa=fx(fitskewdipole);
bya=by(fitskewdipole);
aya=ay(fitskewdipole);
gya=gy(fitskewdipole);
fya=fy(fitskewdipole);






%% Horizontal response matrix derivative respect to quads
sin2k=sin(2*sKlxq)./sKxq;
cos2k=cos(2*sKlxq)-1;

I_k0=(bxq/2+gxq./Kxq/2).*Lq+(bxq/2-gxq./Kxq/2).*sin2k/2+(axq./Kxq/2).*cos2k;
I_ks2=(-cos2k+(axq./bxq).*(sin2k-2*Lq))./Kxq/2;
I_kc2=I_k0+(1./bxq./Kxq).*(sin2k-2*Lq)/2;

% I_k0=I_k0.*Deltas(isquad|isdipole);
% I_ks2=I_ks2.*Deltas(isquad|isdipole);
% I_kc2=I_kc2.*Deltas(isquad|isdipole);

%terms
Amplitude_ij=-sqrt(bxb*bxc)/8/sx/sx2;
ph_ij=fxb*ones(1,nc)-ones(nb,1)*fxc;
ph_ik=fxb*ones(1,nfq)-ones(nb,1)*fxq;
ph_jk=fxc'*ones(1,nfq)-ones(nc,1)*fxq;

C_ij1=cos(abs(ph_ij)-pi*Qx);
S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qx);%

C_ik2=cos(2*abs(ph_ik)-2*pi*Qx);
S_ik2=signtilde(ph_ik).*sin(2*abs(ph_ik)-2*pi*Qx);

C_jk2=cos(2*abs(ph_jk)-2*pi*Qx);
S_jk2=signtilde(ph_jk).*sin(2*abs(ph_jk)-2*pi*Qx);


Sigma_ik2=(ones(nb,1)*I_kc2).*S_ik2-(ones(nb,1)*I_ks2).*C_ik2;
Gamma_ik2=(ones(nb,1)*I_ks2).*S_ik2+(ones(nb,1)*I_kc2).*C_ik2;
Sigma_jk2=(ones(nc,1)*I_kc2).*S_jk2-(ones(nc,1)*I_ks2).*C_jk2;
Gamma_jk2=(ones(nc,1)*I_ks2).*S_jk2+(ones(nc,1)*I_kc2).*C_jk2;

% the complete thing
C=Amplitude_ij.*C_ij1;
S=Amplitude_ij.*S_ij1;

termC1=repmat(permute(Gamma_ik2,[1 3 2]),[1 nc 1]);
termC2=repmat(permute(Gamma_jk2,[3 1 2]),[nb 1 1]);
termC3=repmat(permute(2*I_k0*cx^2,[1 3 2]),[nb nc 1]);
termC=repmat(C,[1 1 nfq]).*(termC1+termC2+termC3);

termS1=repmat(permute(Sigma_ik2+2*sx2*(ones(nb,1)*I_k0).*heaviside(ph_ik),[1 3 2]),[1 nc 1]);
termS2=repmat(permute(Sigma_jk2+2*sx2*(ones(nc,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
termS3=repmat(permute(sx2*I_k0,[1 3 2]),[nb nc 1]).*repmat(-signtilde(ph_ij),[1 1 nfq]);
termS=repmat(S,[1 1 nfq]).*(termS1-termS2+termS3);

dX_dq=termC+termS;

%% Horizontal dispersion derivative respect to quads
sink=sin(sKdlx)./sKdx;
cosk=cos(sKdlx)-1;
I_js1=-cosk./Kdx./bxd;
I_jc1=sink+cosk.*axd./Kdx./bxd;
dI_js1=(sink.*Ld/2+cosk./Kdx)./Kdx./bxd;
dI_jc1=Ld.*(1+cosk)./Kdx/2-sink./Kdx/2-axd.*dI_js1;

% dI_js1=delpa_pdI_js1.*Deltas(isdipole);
% dI_jc1=dI_jc1.*Deltas(isdipole);

h_j=angd./Ld;
%terms
Amplitude_ij=-sqrt(bxb*bxd)/8/sx/sx2;
ph_ij=fxb*ones(1,ndq)-ones(nb,1)*fxd;
ph_ik=fxb*ones(1,nfq)-ones(nb,1)*fxq;
ph_jk=fxd'*ones(1,nfq)-ones(ndq,1)*fxq;

C_ij1=cos(abs(ph_ij)-pi*Qx);
S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qx);%

C_ik2=cos(2*abs(ph_ik)-2*pi*Qx);
S_ik2=signtilde(ph_ik).*sin(2*abs(ph_ik)-2*pi*Qx);

C_jk2=cos(2*abs(ph_jk)-2*pi*Qx);
S_jk2=signtilde(ph_jk).*sin(2*abs(ph_jk)-2*pi*Qx);


Sigma_ik2=(ones(nb,1)*I_kc2).*S_ik2-(ones(nb,1)*I_ks2).*C_ik2;
Gamma_ik2=(ones(nb,1)*I_ks2).*S_ik2+(ones(nb,1)*I_kc2).*C_ik2;
Sigma_jk2=(ones(ndq,1)*I_kc2).*S_jk2-(ones(ndq,1)*I_ks2).*C_jk2;
Gamma_jk2=(ones(ndq,1)*I_ks2).*S_jk2+(ones(ndq,1)*I_kc2).*C_jk2;

% the complete thing R
C=Amplitude_ij.*C_ij1.*(ones(nb,1)*(I_jc1.*h_j));
S=Amplitude_ij.*S_ij1.*(ones(nb,1)*(I_jc1.*h_j));

termC1=repmat(permute(Gamma_ik2,[1 3 2]),[1 ndq 1]);
termC2=repmat(permute(Gamma_jk2,[3 1 2]),[nb 1 1]);
termC3=repmat(permute(2*I_k0*cx^2,[1 3 2]),[nb ndq 1]);
termRC=repmat(C,[1 1 nfq]).*(termC1+termC2+termC3);

termS1=repmat(permute(Sigma_ik2+2*sx2*(ones(nb,1)*I_k0).*heaviside(ph_ik),[1 3 2]),[1 ndq 1]);
termS2=repmat(permute(Sigma_jk2+2*sx2*(ones(ndq,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
termS3=repmat(permute(sx2*I_k0,[1 3 2]),[nb ndq 1]).*repmat(signtilde(ph_ij),[1 1 nfq]);
termRS=repmat(S,[1 1 nfq]).*(termS1-termS2-termS3);


% the complete thing T
C=Amplitude_ij.*S_ij1.*(ones(nb,1)*(I_js1.*h_j));
S=Amplitude_ij.*C_ij1.*(ones(nb,1)*(I_js1.*h_j));

termTC=repmat(C,[1 1 nfq]).*(termC1-termC2+termC3);

termS2=repmat(permute(-Sigma_jk2+2*sx2*(ones(ndq,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
termTS=repmat(S,[1 1 nfq]).*(termS2-termS1+termS3);

%%% below this point test replace ndq with nbd to solve size errors

% %terms
% Amplitude_ij=-sqrt(bxb*bxd)/8/sx/sx2;
% ph_ij=fxb*ones(1,nbd)-ones(nb,1)*fxd;
% ph_ik=fxb*ones(1,nfq)-ones(nb,1)*fxq;
% ph_jk=fxd'*ones(1,nfq)-ones(nbd,1)*fxq;
% 
% C_ij1=cos(abs(ph_ij)-pi*Qx);
% S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qx);%
% 
% C_ik2=cos(2*abs(ph_ik)-2*pi*Qx);
% S_ik2=signtilde(ph_ik).*sin(2*abs(ph_ik)-2*pi*Qx);
% 
% C_jk2=cos(2*abs(ph_jk)-2*pi*Qx);
% S_jk2=signtilde(ph_jk).*sin(2*abs(ph_jk)-2*pi*Qx);
% 
% 
% Sigma_ik2=(ones(nb,1)*I_kc2).*S_ik2-(ones(nb,1)*I_ks2).*C_ik2;
% Gamma_ik2=(ones(nb,1)*I_ks2).*S_ik2+(ones(nb,1)*I_kc2).*C_ik2;
% Sigma_jk2=(ones(nbd,1)*I_kc2).*S_jk2-(ones(nbd,1)*I_ks2).*C_jk2;
% Gamma_jk2=(ones(nbd,1)*I_ks2).*S_jk2+(ones(nbd,1)*I_kc2).*C_jk2;
% 
% % the complete thing R
% C=Amplitude_ij.*C_ij1.*(ones(nb,1)*(I_jc1.*h_j));
% S=Amplitude_ij.*S_ij1.*(ones(nb,1)*(I_jc1.*h_j));
% 
% termC1=repmat(permute(Gamma_ik2,[1 3 2]),[1 nbd 1]);
% termC2=repmat(permute(Gamma_jk2,[3 1 2]),[nb 1 1]);
% termC3=repmat(permute(2*I_k0*cx^2,[1 3 2]),[nb nbd 1]);
% termRC=repmat(C,[1 1 nfq]).*(termC1+termC2+termC3);
% 
% termS1=repmat(permute(Sigma_ik2+2*sx2*(ones(nb,1)*I_k0).*heaviside(ph_ik),[1 3 2]),[1 nbd 1]);
% termS2=repmat(permute(Sigma_jk2+2*sx2*(ones(nbd,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
% termS3=repmat(permute(sx2*I_k0,[1 3 2]),[nb nbd 1]).*repmat(signtilde(ph_ij),[1 1 nfq]);
% termRS=repmat(S,[1 1 nfq]).*(termS1-termS2-termS3);
% 
% 
% 
% % the complete thing T
% C=Amplitude_ij.*S_ij1.*(ones(nb,1)*(I_js1.*h_j));
% S=Amplitude_ij.*C_ij1.*(ones(nb,1)*(I_js1.*h_j));
% 
% termTC=repmat(C,[1 1 nfq]).*(termC1-termC2+termC3);
% 
% termS2=repmat(permute(-Sigma_jk2+2*sx2*(ones(nbd,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
% termTS=repmat(S,[1 1 nfq]).*(termS2-termS1+termS3);

% the dipoles term
Amplitude_ik=sqrt(bxb*(bxd))/2/sx;
if nbd>0 & nbd<nfq
%    warning('CalcRespXXRespMat_thick_V2: problem if more dipoles than quadrupoles!')
    dipoleterm=[zeros(nb,nfq-nbd) (ones(nb,1)*h_j).*Amplitude_ik.*((ones(nb,1)*dI_js1).*S_ij1+(ones(nb,1)*dI_jc1).*C_ij1)];
else
    dipoleterm=[zeros(nb,nfq)];
end

dDx_dq=squeeze(sum(termRC+termRS+termTC+termTS,2))+dipoleterm;

term_direct=(ones(nb,1)*I_jc1).*C_ij1+(ones(nb,1)*I_js1).*S_ij1;
term_alpha=-(ones(nb,1)*(h_j.*cosk./Kdx.*tan(e1))).*C_ij1;
term_rho2=2*((ones(nb,1)*dI_js1).*S_ij1+(ones(nb,1)*dI_jc1).*C_ij1).*(ones(nb,1)*h_j).^2;
dDx_db=Amplitude_ik.*(term_direct+term_alpha+term_rho2);

Dx=sum(Amplitude_ik.*term_direct.*(ones(nb,1)*h_j),2);

%% Vertical response matrix derivative
sin2k=sin(2*sKlyq)./sKyq;
cos2k=cos(2*sKlyq)-1;


I_k0=(byq/2+gyq./Kyq/2).*Lq+(byq/2-gyq./Kyq/2).*sin2k/2+(ayq./Kyq/2).*cos2k;
I_ks2=(-cos2k+(ayq./byq).*(sin2k-2*Lq))./Kyq/2;
I_kc2=I_k0+(1./byq./Kyq).*(sin2k-2*Lq)/2;

% I_k0=I_k0.*Deltas(isquad|isdipole);
% I_ks2=I_ks2.*Deltas(isquad|isdipole);
% I_kc2=I_kc2.*Deltas(isquad|isdipole);

%terms
Amplitude_ij=sqrt(byb*byc)/8/sy/sy2;
ph_ij=fyb*ones(1,nc)-ones(nb,1)*fyc;
ph_ik=fyb*ones(1,nfq)-ones(nb,1)*fyq;
ph_jk=fyc'*ones(1,nfq)-ones(nc,1)*fyq;

C_ij1=cos(abs(ph_ij)-pi*Qy);
S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qy);%

C_ik2=cos(2*abs(ph_ik)-2*pi*Qy);
S_ik2=signtilde(ph_ik).*sin(2*abs(ph_ik)-2*pi*Qy);

C_jk2=cos(2*abs(ph_jk)-2*pi*Qy);
S_jk2=signtilde(ph_jk).*sin(2*abs(ph_jk)-2*pi*Qy);


Sigma_ik2=(ones(nb,1)*I_kc2).*S_ik2-(ones(nb,1)*I_ks2).*C_ik2;
Gamma_ik2=(ones(nb,1)*I_ks2).*S_ik2+(ones(nb,1)*I_kc2).*C_ik2;
Sigma_jk2=(ones(nc,1)*I_kc2).*S_jk2-(ones(nc,1)*I_ks2).*C_jk2;
Gamma_jk2=(ones(nc,1)*I_ks2).*S_jk2+(ones(nc,1)*I_kc2).*C_jk2;

% the complete thing
C=Amplitude_ij.*C_ij1;
S=Amplitude_ij.*S_ij1;

termC1=repmat(permute(Gamma_ik2,[1 3 2]),[1 nc 1]);
termC2=repmat(permute(Gamma_jk2,[3 1 2]),[nb 1 1]);
termC3=repmat(permute(2*I_k0*cy^2,[1 3 2]),[nb nc 1]);
termC=repmat(C,[1 1 nfq]).*(termC1+termC2+termC3);

termS1=repmat(permute(Sigma_ik2+2*sy2*(ones(nb,1)*I_k0).*heaviside(ph_ik),[1 3 2]),[1 nc 1]);
termS2=repmat(permute(Sigma_jk2+2*sy2*(ones(nc,1)*I_k0).*heaviside(ph_jk),[3 1 2]),[nb 1 1]);
termS3=repmat(permute(sy2*I_k0,[1 3 2]),[nb nc 1]).*repmat(-signtilde(ph_ij),[1 1 nfq]);
termS=repmat(S,[1 1 nfq]).*(termS1-termS2+termS3);

dY_dq=termC+termS;

%% Sxy and Sxy skew term
spos=sin(pi*(Qx+Qy));
sneg=sin(pi*(Qx-Qy));

Order_mj=heaviside(repmat(bpmlist,[1,nc,nfsq])-repmat(reshape(fitskewquadlist,[1,1,nfsq]),[nb,nc,1]));
Order_mw=heaviside(repmat(cmlist,[nb,1,nfsq])-repmat(reshape(fitskewquadlist,[1,1,nfsq]),[nb,nc,1]));
Order_jw=heaviside(repmat(cmlist,[nb,1,nfsq])-repmat(bpmlist,[1,nc,nfsq]));

dfx_wj=repmat(fxb,[1,nc,nfsq])-repmat(fxc,[nb,1,nfsq])+2*pi*Qx*(Order_jw);
dfx_mj=repmat(fxb,[1,nc,nfsq])-repmat(reshape(fxs,[1,1,nfsq]),[nb,nc,1])+2*pi*Qx*not(Order_mj);
dfx_mw=repmat(fxc,[nb,1,nfsq])-repmat(reshape(fxs,[1,1,nfsq]),[nb,nc,1])+2*pi*Qx*not(Order_mw);

tfx_wj=dfx_wj-pi*Qx;
tfx_mj=dfx_mj-pi*Qx;
tfx_mw=dfx_mw-pi*Qx;

dfy_wj=repmat(fyb,[1,nc,nfsq])-repmat(fyc,[nb,1,nfsq])+2*pi*Qy*(Order_jw);
dfy_mj=repmat(fyb,[1,nc,nfsq])-repmat(reshape(fys,[1,1,nfsq]),[nb,nc,1])+2*pi*Qy*not(Order_mj);
dfy_mw=repmat(fyc,[nb,1,nfsq])-repmat(reshape(fys,[1,1,nfsq]),[nb,nc,1])+2*pi*Qy*not(Order_mw);

tfy_wj=dfy_wj-pi*Qy;
tfy_mj=dfy_mj-pi*Qy;
tfy_mw=dfy_mw-pi*Qy;

bys_re=repmat(reshape(bys.*Ls,[1,1,nfsq]),[nb,nc,1]);%.*Deltas(isskew)
bxs_re=repmat(reshape(bxs.*Ls,[1,1,nfsq]),[nb,nc,1]);%.*Deltas(isskew)
Order=Order_mj-Order_mw+Order_jw;


Bxy=sqrt(repmat(bxb*byc,[1,1,nfsq]).*bxs_re.*bys_re);
Byx=sqrt(repmat(byb*bxc,[1,1,nfsq]).*bxs_re.*bys_re);
Aneg=(-cos(tfx_mj-tfy_mj-tfx_wj)/sx+cos(tfx_mw-tfy_mw-tfy_wj)/sy)/sneg;
Apos=(cos(tfx_mj+tfy_mj-tfx_wj)/sx+cos(tfx_mw+tfy_mw+tfy_wj)/sy)/spos;
dYX_ds=Byx/8.*(Apos+Aneg);

Aneg=(cos(tfx_mj-tfy_mj+tfy_wj)/sy-cos(tfx_mw-tfy_mw+tfx_wj)/sx)/sneg;
Apos=(cos(tfx_mj+tfy_mj-tfy_wj)/sy+cos(tfx_mw+tfy_mw+tfx_wj)/sx)/spos;
dXY_ds=Bxy/8.*(Apos+Aneg);


%XY=sum(dXY_ds.*repmat(reshape(Ks,[1,1,nfsq]),[nb,nc,1]),3);
%YX=sum(dYX_ds.*repmat(reshape(Ks,[1,1,nfsq]),[nb,nc,1]),3);

%% Dy and Dy dipole and skew dipole term


%sink=sin(sKlxfb)./sKxfb;
%cosk=cos(sKlxfb)-1;
%I_js1=-cosk./Kxfb./bxfb;
%I_jc1=sink+cosk.*axfb./Kxfb./bxfb;
I_jc1=cos(angd/2)- axd0.*(Ld/2).*sin(angd/2)/2./(angd/2)./bxfb;
I_js1=(Ld/2).*sin(angd/2)/2./(angd/2)./bxfb;

Amplitude_ik=sqrt(bxb*(bxfb))/2/sx;
ph_ij=fxb*ones(1,nfd)-ones(nb,1)*fxfb;
C_ij1=cos(abs(ph_ij)-pi*Qx);
S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qx);
term_direct=(ones(nb,1)*I_jc1).*C_ij1+(ones(nb,1)*I_js1).*S_ij1;
dDx_db=Amplitude_ik.*term_direct;

% senza integrali
% AmplitudeCenter_ik=sqrt(bxb*(bxfbcenter))/2/sx;
% ph_ijcenter=fxbcenter*ones(1,nfd)-ones(nb,1)*fxfbcenter;
% C_ij1center=cos(abs(ph_ijcenter)-pi*Qx);
% dDx_db=AmplitudeCenter_ik.*C_ij1center;

Dx=sum(dDx_db.*(ones(nb,1)*(angd./2)),2);

sink=sin(sKlya)./sKya;
cosk=cos(sKlya)-1;
I_js1=-cosk./Kya./bya;
I_jc1=sink+cosk.*ayfb./Kya./bya;
Amplitude_ik=sqrt(byb*(bya))/2/sy;
ph_ij=fyb*ones(1,nfsd)-ones(nb,1)*fya;
C_ij1=cos(abs(ph_ij)-pi*Qy);
S_ij1=signtilde(ph_ij).*sin(abs(ph_ij)-pi*Qy);
term_direct=(ones(nb,1)*I_jc1).*C_ij1+(ones(nb,1)*I_js1).*S_ij1;
dDy_da=Amplitude_ik.*term_direct;

% empty skew dipole Parmas, size errors, set to zero. This is exactly what
% I would like for the dispersion fit!
%dDy_da=zeros(size(dDx_db));


%% Dy and Dy skew quad term
tfx_mjsq=squeeze(tfx_mj(:,1,:));
tfy_mjsq=squeeze(tfy_mj(:,1,:));

ssh=sin(sKlxs).*sinh(sKlxs);
cch=cos(sKlxs).*cosh(sKlxs);
sch=sin(sKlxs).*cosh(sKlxs);
csh=cos(sKlxs).*sinh(sKlxs);


TSx=dys./sqrt(bxs)./Kxs/2.*(ssh-cch+1)+dpys./sqrt(bxs)./sKlxs./Kxs/2.*(sch-csh);

TCx=dys.*sqrt(bxs)./sKxs/2.*(sch+csh)-axs.*dpys./sqrt(bxs)./sKlxs./Kxs/2.*(sch-csh)-axs.*dys./sqrt(bxs)./Kxs/2.*(ssh-cch+1)+dpys.*sqrt(bxs)./Kxs/2.*(ssh+cch-1);
TSy=dxs./sqrt(bys)./Kxs/2.*(ssh+cch-1)+dpxs./sqrt(bys)./sKlxs./Kxs/2.*(sch-csh);
TCy=dxs.*sqrt(bys)./sKxs/2.*(sch+csh)-ays.*dpxs./sqrt(bys)./sKlxs./Kxs/2.*(sch-csh)-ays.*dxs./sqrt(bys)./Kxs/2.*(ssh+cch-1)+dpxs.*sqrt(bys)./Kxs/2.*(ssh-cch+1);

dDx_ds=(sqrt(bxb)*TCx.*cos(tfx_mjsq)+sqrt(bxb)*TSx.*sin(tfy_mjsq))/2/sx;
dDy_ds=(sqrt(byb)*TCy.*cos(tfy_mjsq)+sqrt(byb)*TSy.*sin(tfx_mjsq))/2/sy;

% thin case pluse second order term
%dDx_ds=sqrt(bxb*bxs).*repmat((dys+Ks.*Ls.*bys.*dxs*cy/2/sy).*Ls,[nb,1]).*cos(tfx_mjsq)/2/sx;%.*Deltas(isskew)
%dDy_ds=sqrt(byb*bys).*repmat((dxs+Ks.*Ls.*bxs.*dys*cx/2/sx).*Ls,[nb,1]).*cos(tfy_mjsq)/2/sy;%.*Deltas(isskew)
end

function y=heaviside(x)
y=sign(x).*(1+signtilde(x))/2;
end
function y=signtilde(x)
y=sign(x)-double(x==0);
end



