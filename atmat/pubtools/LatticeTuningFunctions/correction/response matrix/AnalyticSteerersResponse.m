function [rmh,rmv]=AnalyticSteerersResponse(mach,bpmidx,coridx)
%SEMRDT compute resonance driving terms at BPM locations
%
% [rmh,rmv]=AnalyticSteerersResponse(mach,bpmidx,coridx)%
% mach    : AT lattice
% bpmindx : BPM indexes
% coridx : correctors indexes
% 
% rmh    : rm D Hbpm /H Vcor (=1)
% rmv    : rm D Vbpm /D Vcor (=1)
%
% to obtain orbit response for a given set of coridx strengths (KL)
% 
% x=rmh.*k1n.*Lcor
% y=rmv.*k1s.*Lcor
% 
% this function is inspired by semrdtresp by L.Farvacque
% 
%see also: atavedata_mod

nb=length(bpmidx);
ns=length(coridx);

% Compute optics

[refpts,ii,kl]=unique([coridx bpmidx length(mach)+1]);
jsk=kl(1:ns);
jbpm=kl(ns+(1:nb));
jend=kl(end);
[vdata,avebeta,avemu,avedisp]=atavedata_mod(mach,0,refpts);
mtunes=vdata(jend).mu;
if ~isempty(find(avebeta<0))
    bx=arrayfun(@(a)a.beta(1),vdata);
    by=arrayfun(@(a)a.beta(2),vdata);
    avebeta=[bx,by];
    dx=arrayfun(@(a)a.Dispersion(1),vdata);
    dy=arrayfun(@(a)a.Dispersion(3),vdata);
    avedisp=[dx,dy];
    warning('on','all');
    warning('negative data in AVEBETA! using beta at entrance!')
    save('failingavebetalattice.mat','mach','bpmidx','skewidx')
    warning('off','all');
end
% Extract parameters

bpm.beta=cat(1,vdata(jbpm).beta);
bpm.phase=cat(1,vdata(jbpm).mu);
bpm.disp=bpm.phase;
bpm.disp(:,1)=arrayfun(@(a)a.Dispersion(1),vdata(jbpm));
bpm.disp(:,2)=arrayfun(@(a)a.Dispersion(3),vdata(jbpm));

cor.beta=avebeta(jsk,:);
cor.phase=avemu(jsk,:);
cor.disp=avedisp(jsk,:);

% cor.beta=cat(1,vdata(jsk).beta);
% cor.phase=cat(1,vdata(jsk).mu);
% cor.disp=cor.phase;
% cor.disp(:,1)=arrayfun(@(a)a.Dispersion(1),vdata(jsk));
% cor.disp(:,2)=arrayfun(@(a)a.Dispersion(3),vdata(jsk));


% Compute terms rm dx/dhcor
Lh=atgetfieldvalues(mach,coridx,'Length');   
Lh(Lh==0)=1;% thin lens magnets
lengthsHmat=repmat(Lh',length(bpmidx),1);

jsqbx=real(sqrt(repmat(bpm.beta(:,1),1,length(cor.beta(:,1))).*repmat(cor.beta(:,1)',length(bpm.beta(:,1)),1)));
dphix=(repmat(bpm.phase(:,1),1,length(cor.phase(:,1)))-repmat(cor.phase(:,1)',length(bpm.phase(:,1)),1));
neg=(dphix < 0);
dphix(neg)=dphix(neg)+mtunes(1);
dphix=abs(dphix);%-mtunes(1);
dis=(repmat(bpm.disp(:,1),1,length(cor.disp(:,1)))-repmat(cor.disp(:,1)',length(bpm.disp(:,1)),1))/findspos(mach,length(mach)+1)/mcf(mach);

re=jsqbx.*cos(dphix);
denom=2*sin(mtunes(1));

rmh=(re/denom+dis).*lengthsHmat;

% Compute terms rm dy/dvcor
Lv=atgetfieldvalues(mach,coridx,'Length');   
Lv(Lv==0)=1;% thin lens magnets
lengthsVmat=repmat(Lv',length(bpmidx),1);
jsqbz=real(sqrt(repmat(bpm.beta(:,2),1,length(cor.beta(:,2))).*repmat(cor.beta(:,2)',length(bpm.beta(:,2)),1)));
dphiz=(repmat(bpm.phase(:,2),1,length(cor.phase(:,2)))-repmat(cor.phase(:,2)',length(bpm.phase(:,2)),1));
neg=(dphiz < 0);
dphiz(neg)=dphiz(neg)+mtunes(2);
dphiz=abs(dphiz);%-mtunes(2);
dis=(repmat(bpm.disp(:,2),1,length(cor.disp(:,2)))-repmat(cor.disp(:,2)',length(bpm.disp(:,2)),1))/findspos(mach,length(mach)+1)/mcf(mach);

re=jsqbz.*cos(dphiz);
denom=2*sin(mtunes(2));

rmv=(re/denom+dis).*lengthsVmat;


end
