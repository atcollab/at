function [d0,de,dc]=DisplayCorrectionEffect(...
    r0,...
    rerr,...
    rcor,...
    inCODe,...
    inCODc,...
    refpts,...
    indHCor,...
    indVCor,...
    indQCor,...
    indSCor)
% [d0,de,dc]=DisplayCorrectionEffect(...
%     r0,...        1) reference lattice
%     rerr,...      2) lattice with errors
%     rcor,...      3) corrected lattice
%     inCODe,...    4) initial COD guess for rerr
%     inCODc,...    5) initial COD guess for rcor
%     refpts,...    6) reference points for computation
%     indHCor,...   7) hor. steerers indexes
%     indVCor,...   8) ver. steerers indexes 
%     indQCor,...   9) normal quad. correctors indexes
%     indSCor)     10) skew quad. correctors indexes
%
% dispaly correction effect.
%

compute_emittances=true;

disp(' --- model lattice data --- ')

d0=getdatalattice(r0,inCODe*0,refpts,indHCor,indVCor,indQCor,indSCor,compute_emittances);

disp(' --- errors lattice data --- ')
de=getdatalattice(rerr,inCODe,refpts,indHCor,indVCor,indQCor,indSCor,compute_emittances);

disp(' --- corrected lattice data --- ')
dc=getdatalattice(rcor,inCODc,refpts,indHCor,indVCor,indQCor,indSCor,compute_emittances);


% print out correction effect
oudataforma='%3.2e';
oudataformatune='%2.3f';
oudataformaemit='%3.3f';
oudataformabeta='%2.1f';

%
disp('--------------------------------------------------');
disp('------ total std corrector values applyed --------')
disp('------                                    --------');
disp(['       HK (' num2str(length(d0.ch)) ') [1/m]: ' num2str(std(de.ch),oudataforma) ' -> '  num2str(std(dc.ch),oudataforma) ]);
disp(['       VK (' num2str(length(d0.cv)) ') [1/m]: ' num2str(std(de.cv),oudataforma) ' -> '  num2str(std(dc.cv),oudataforma) ]);
disp(['       SK (' num2str(length(d0.cs)) ') [1/m2]: ' num2str(std(de.cs),oudataforma) ' -> '  num2str(std(dc.cs),oudataforma) ]);
disp(['       QK (' num2str(length(d0.cq)) ') [1/m2]: ' num2str(std(de.cq-d0.cq),oudataforma) ' -> '  num2str(std(dc.cq-d0.cq),oudataforma) ]);
if nargin==8
    ch=dc.ch.*d0.Lh;
    cv=dc.cv.*d0.Lv;
    cq=dc.cq.*d0.Lq;
    cs=dc.cs.*d0.Ls;
    che=de.ch.*d0.Lh;
    cve=de.cv.*d0.Lv;
    cqe=de.cq.*d0.Lq;
    cse=de.cs.*d0.Ls;
    disp(['       HKL (' num2str(length(ch)) ') [rad]: ' num2str(std(che),oudataforma) ' -> '  num2str(std(ch),oudataforma) ]);
    disp(['       VKL (' num2str(length(cv)) ') [rad]: ' num2str(std(cve),oudataforma) ' -> '  num2str(std(cv),oudataforma) ]);
    disp(['       SKL (' num2str(length(cs)) ') [1/m]: ' num2str(std(cse),oudataforma) ' -> '  num2str(std(cs),oudataforma) ]);
    disp(['       QKL (' num2str(length(cq)) ') [1/m]: ' num2str(std(cqe-d0.cq),oudataforma) ' -> '  num2str(std(cq-d0.cq),oudataforma) ]);
    
    Brho=getBrho(r0);
    
    ch=dc.ch.*d0.Lh.*Brho;
    cv=dc.cv.*d0.Lv.*Brho;
    cq=dc.cq.*d0.Lq.*Brho;
    cs=dc.cs.*d0.Ls.*Brho;
    che=de.ch.*d0.Lh.*Brho;
    cve=de.cv.*d0.Lv.*Brho;
    cqe=de.cq.*d0.Lq.*Brho;
    cse=de.cs.*d0.Ls.*Brho;
    disp(['       HKLBrho (' num2str(length(ch)) ') [Tm]: ' num2str(std(che),oudataforma) ' -> '  num2str(std(ch),oudataforma) ]);
    disp(['       VKLBrho (' num2str(length(cv)) ') [Tm]: ' num2str(std(cve),oudataforma) ' -> '  num2str(std(cv),oudataforma) ]);
    disp(['       SKLBrho (' num2str(length(cs)) ') [T]: ' num2str(std(cse),oudataforma) ' -> '  num2str(std(cs),oudataforma) ]);
    disp(['       QKLBrho (' num2str(length(cq)) ') [T]: ' num2str(std(cqe-d0.cq),oudataforma) ' -> '  num2str(std(cq-d0.cq),oudataforma) ]);
end
disp('------                                    --------');
disp('------    residual orbit and dispersion   --------')
disp('------                                    --------');
disp(['       OH (' num2str(length(d0.monh)) ') [m]: ' num2str(std(de.monh),oudataforma) ' -> '  num2str(std(dc.monh),oudataforma) ]);
disp(['       OV (' num2str(length(d0.monv)) ') [m]: ' num2str(std(de.monv),oudataforma) ' -> '  num2str(std(dc.monv),oudataforma) ]);
disp(['       DH (' num2str(length(d0.dish)) ') [m]: ' num2str(std(de.dish-d0.dish),oudataforma) ' -> '  num2str(std(dc.dish-d0.dish),oudataforma) ]);
disp(['       DV (' num2str(length(d0.disv)) ') [m]:' num2str(std(de.disv),oudataforma) ' -> '  num2str(std(dc.disv),oudataforma) ]);
disp(['       BBH (' num2str(length(d0.bbh)) ') %: ' num2str(std((de.bbh-d0.bbh)./d0.bbh)*100,oudataformabeta) ' -> '  num2str(std((dc.bbh-d0.bbh)./d0.bbh)*100,oudataformabeta) ]);
disp(['       BBV (' num2str(length(d0.bbv)) ') %: ' num2str(std((de.bbv-d0.bbv)./d0.bbv)*100,oudataformabeta) ' -> '  num2str(std((dc.bbv-d0.bbv)./d0.bbv)*100,oudataformabeta) ]);
disp(['       PhH (' num2str(length(d0.mh)) ') : ' num2str(std((de.mh-d0.mh)),oudataforma) ' -> '  num2str(std((dc.mh-d0.mh)),oudataforma) ]);
disp(['       PhV (' num2str(length(d0.mv)) ') : ' num2str(std((de.mv-d0.mv)),oudataforma) ' -> '  num2str(std((dc.mv-d0.mv)),oudataforma) ]);
disp('------                                    --------');
disp('------        tune and emittance          --------')
disp('------                                    --------');
disp(['       Qx [' num2str(d0.tune(1),oudataformatune) ']: ' num2str(de.tune(1),oudataformatune) ' -> '  num2str(dc.tune(1),oudataformatune) ]);
disp(['       Qy [' num2str(d0.tune(2),oudataformatune) ']: ' num2str(de.tune(2),oudataformatune) ' -> '  num2str(dc.tune(2),oudataformatune) ]);
disp(['       Cx [' num2str(d0.crom(1),oudataformatune) ']: ' num2str(de.crom(1),oudataformatune) ' -> '  num2str(dc.crom(1),oudataformatune) ]);
disp(['       Cy [' num2str(d0.crom(2),oudataformatune) ']: ' num2str(de.crom(2),oudataformatune) ' -> '  num2str(dc.crom(2),oudataformatune) ]);
if compute_emittances
    disp(['       EX [' num2str(d0.modemittance(1)*1e12,oudataformaemit) ' pm]: ' num2str(de.modemittance(1)*1e12,oudataformaemit) ' -> '  num2str(dc.modemittance(1)*1e12,oudataformaemit) ]);
    disp(['       EY [' num2str(d0.modemittance(2)*1e12,oudataformaemit) 'pm]: ' num2str(de.modemittance(2)*1e12,oudataformaemit) ' -> '  num2str(dc.modemittance(2)*1e12,oudataformaemit) ]);
end
disp('------                                    --------');
disp('--------------------------------------------------');


return


function a=getdatalattice(r0,inCOD,refpts,indHCor,indVCor,indCorQuads,indSCor,emitok)

warning('off','all'); % mcf,atx, findorbit6,... generates warnings

alpha=mcf(r0);
indrfc=find(atgetcells(r0,'Frequency'));
            
% get initial orbit
o=findorbit6Err(r0,refpts,inCOD);
a.monh=o(1,:);
a.monv=o(3,:);
d=finddispersion6Err(r0,refpts,indrfc,alpha,1e-4,inCOD);
a.dish=d(1,:);
a.disv=d(3,:);

if emitok
    try
        [~,b0]=atx(r0,0,1:length(r0));
    catch exc
        getReport(exc,'extended');
        warning('atx failed');
        b0.modemittance=[NaN NaN];
        b0.fulltunes=[NaN NaN];
    end
    
    a.tune=b0.fulltunes;
    a.modemittance= b0.modemittance;
end

[l,t,a.crom]=atlinopt(r0,0,refpts);

if ~emitok
    a.tune=t;
end

a.bbh=arrayfun(@(s)s.beta(1),l);
a.bbv=arrayfun(@(s)s.beta(2),l);
a.mh=arrayfun(@(s)s.mu(1),l);
a.mv=arrayfun(@(s)s.mu(2),l);

a.Lh=getcellstruct(r0,'Length',indHCor);
a.Lv=getcellstruct(r0,'Length',indVCor);
a.Lq=getcellstruct(r0,'Length',indCorQuads);
a.Ls=getcellstruct(r0,'Length',indSCor);
a.ch=getcellstruct(r0,'PolynomB',indHCor,1,1);
a.cv=getcellstruct(r0,'PolynomA',indVCor,1,1);
a.cq=getcellstruct(r0,'PolynomB',indCorQuads,1,2);
a.cs=getcellstruct(r0,'PolynomA',indSCor,1,2);

warning('on','all');

return
