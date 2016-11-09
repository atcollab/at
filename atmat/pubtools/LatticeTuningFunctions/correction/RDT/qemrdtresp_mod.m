function [fx,fz,qcor]=qemrdtresp_mod(mach,bpmidx,qcoridx)
%QEMRDTRESP  compute resonance driving terms at BPM locations
%
%[fx,fz,qcor]=semrdtresp_mod(mach,bpmidx,qcoridx)
%
% mach    : AT lattice
% bpmindx : BPM indexes
% qcoridx : skew quadrupole indexes
% 
% fx    : f2000 RDT 
% fx    : f0020 RDT
% qcor  : qcor.beta, qcor.phase beta and phase at the quad index (averaged) 
% 
% to obtain rdt for a given set of skewidx strengths (KL)
% 
% f1001=f1.*k1.*Lquad
% f1010=f2.*k1.*Lquad
% 
% this function is an exact copy of semrdtresp by L.Farvacque
% 
%see also: atavedata_mod

nbpm=length(bpmidx);
nqcor=length(qcoridx);

% Compute optics

[refpts,~,kl]=unique([qcoridx bpmidx length(mach)+1]);
jcor=kl(1:nqcor);
jbpm=kl(nqcor+(1:nbpm));
jend=kl(end);
[vdata,avebeta,avemu]=atavedata_mod(mach,0,refpts);
mtx=vdata(jend).mu(1);
mtz=vdata(jend).mu(2);

% Extract parameters
if nargout >= 3
    qcor.beta=avebeta(jcor,:);
    qcor.phase=avemu(jcor,:);
end

% Compute terms

dphix=dphase(avemu(jbpm,1),avemu(jcor,1),mtx);
dphiz=dphase(avemu(jbpm,2),avemu(jcor,2),mtz);

fx=-avebeta(jcor,ones(1,nbpm))'.*complex(cos(2*dphix),sin(2*dphix))./...
    (1-complex(cos(2*mtx),sin(2*mtx)))/8;

fz= avebeta(jcor,2*ones(1,nbpm))'.*complex(cos(2*dphiz),sin(2*dphiz))./...
    (1-complex(cos(2*mtz),sin(2*mtz)))/8;

    function dph=dphase(phib,phik,mtune)
        nb=length(phib);
        nk=length(phik);
        dph=phik(:,ones(nb,1))'-phib(:,ones(1,nk));
        neg=(dph < 0);
        dph(neg)=dph(neg)+mtune;
    end

end
