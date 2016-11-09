function [f1,f2,skew]=semrdtresp_mod(mach,bpmidx,skewidx)
%SEMRDT compute resonance driving terms at BPM locations
%
%[f1,f2,skew]=semrdtresp_mod(mach,bpmidx,skewidx)
%
% mach    : AT lattice
% bpmindx : BPM indexes
% skewidx : skew quadrupole indexes
% 
% f1    : f1001 RDT 
% f2    : f1010 RDT
% skew  : skew.beta skew.phase beta and phase at the skew index (averaged) 
% 
% to obtain rdt for a given set of skewidx strengths (KL)
% 
% f1001=f1.*k1s.*Lskew
% f1010=f2.*k1s.*Lskew
% 
% this function is an exact copy of semrdtresp by L.Farvacque
% 
%see also: atavedata_mod

nb=length(bpmidx);
ns=length(skewidx);

% Compute optics

[refpts,ii,kl]=unique([skewidx bpmidx length(mach)+1]);
jsk=kl(1:ns);
jbpm=kl(ns+(1:nb));
jend=kl(end);
[vdata,avebeta,avemu]=atavedata_mod(mach,0,refpts);
mtunes=vdata(jend).mu;
if ~isempty(find(avebeta<0))
    bx=arrayfun(@(a)a.beta(1),vdata);
    by=arrayfun(@(a)a.beta(2),vdata);
    avebeta=[bx,by];
    warning('on','all');
    warning('negative data in AVEBETA! using beta at entrance!')
    save('failingavebetalattice.mat','mach','bpmidx','skewidx')
    warning('off','all');
end
% Extract parameters

bpm.phase=cat(1,vdata(jbpm).mu);

skew.beta=avebeta(jsk,:);
skew.phase=avemu(jsk,:);

% Compute terms

jsqb=real(sqrt(skew.beta(:,1).*skew.beta(:,2)));
[dphix,dphiz]=dphase(bpm.phase,skew.phase',mtunes);

re1=jsqb(:,ones(1,nb))'.*cos(dphix-dphiz);
im1=jsqb(:,ones(1,nb))'.*sin(dphix-dphiz);
t1=mtunes(1)-mtunes(2);

denom1=4*(1-complex(cos(t1),sin(t1)));
f1=complex(re1,im1)/denom1;

re2=jsqb(:,ones(1,nb))'.*cos(dphix+dphiz);
im2=jsqb(:,ones(1,nb))'.*sin(dphix+dphiz);
t2=mtunes(1)+mtunes(2);
denom2=4*(1-complex(cos(t2),sin(t2)));
f2=complex(re2,im2)/denom2;

end

function [dphix,dphiz]=dphase(phib,phik,mtune)
nb=length(phib);
nk=length(phik);
dphix=phik(  ones(nb,1),:)-phib(:,  ones(1,nk));
neg=(dphix < 0);
dphix(neg)=dphix(neg)+mtune(1);
dphiz=phik(2*ones(nb,1),:)-phib(:,2*ones(1,nk));
neg=(dphiz < 0);
dphiz(neg)=dphiz(neg)+mtune(2);
end
