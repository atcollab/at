function [rcor,inCOD,hs,vs]=atdispersionfreesteering(...
    rerr,...
    indBPM,...
    indHCor,...
    indVCor,...
    inCOD,...
    neigSteerer,...
    correctflags,...
    scalefactor,...
    wdisp,...
    ModelRM,...
    reforbit,...
    refdispersion,...
    steererlimit,...
    printouttext)
% function [...
%    rcor,...           1) corrected lattice
%    inCOD,...          2) initial COD (dpp is stored here)
%    hs,vs...           3) required steerers strengths (total)
%    ]=atdispersionfreesteering(...
%     rerr,...          1) AT lattice to correct
%     indBPM,...        2) Nbx1 bpm indexes
%     indHCor,...       3) Nhx1 hor. cor indexes
%     indVCor,...       4) Nvx1 ver. cor indexes
%     inCOD,...         5) 6x1 initial COD guess
%     neigSteerer,...   6) 2xNiter eigenvectors for correction H and V at
%                          each iteration (default: [Nh/2 Nv/2])
%     correctflags,...  7) correct [dpp mean0](default: [true true])
%     scalefactor,...   8) scale factor to correction (default: 1.0)
%     wdisp,...         9) weight dispersion*wdisp and orbit*(1-wdisp)
%                          (default: 0.7)
%     ModelRM,...       10) ModelRM.Orb(H/V)Cor = 6x1 cell of orbit response mat.
%                           ModelRM.Orb(H/V)DPP = 6x1 array of orbit
%                           ModelRM.Disp(H/V)Cor = 6x1 cell of dispersion response mat.
%                           ModelRM.Disp(H/V)DPP = 6x1 array of dispersion
%                           if [] compute RM (default: [])
%     reforbit,...      11) 2xNbpm reference orbit to correct to (default 0*2xNb)
%     refdispersion,... 12) 2xNbpm reference orbit to correct to (default 0*2xNb)
%     steererlimit      13) 2x1 limit of steerers abs(steerer)<steererlimit
%                           (default: [], no limits)
%     printouttext      14) if 1 or true, display rms orbit
%     )
%
% features impelemented:
% limit correctors strengths
% ddp correction
% sum of steerers = 0
% 6D orbit with BPM errors
% initial coordinate
% correction to reference orbit refx refy
% use atsetfieldvalues, atgetcells
%
%
%http://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.3.121001
%
%see also: qemsvd_mod findorbit6Err getresponsematrices



% response matrix kicks
kval=1e-5;
delta=1e-3;

% default arguments
if nargin<14
    printouttext=true;
end
if nargin<13
    steererlimit=[];
end

if nargin<4
    if printouttext
        disp('get BPM and Correctors indexes'); end;
    indBPM=finc(atgetcells(rerr,'Class','Monitor'));
    indHCor=finc(atgetcells(rerr,'iscorH','H'));
    indVCor=finc(atgetcells(rerr,'iscorV','V'));
end

if nargin<5
    inCOD=[0 0 0 0 0 0]';
end

if nargin<6
    neigSteerer=[length(indHCor) length(indVCor)]/2;
end

if nargin<7
    correctflags=[true true];
end

if nargin<8
    if printouttext
        disp(' --- scale set to 1.0'); end;
    scalefactor=1.0;
end

if nargin<9
    if printouttext, disp(' --- alpha=0.7'); end;
    wdisp=0.8;
end

if nargin<10
    if printouttext, disp(' --- computing orbit Response matrix'); end;
    ModelRM=[];
end

if nargin<11
    if printouttext, disp(' --- reference orbit = 0'); end;
    reforbit=zeros(size(indBPM),2);
end

if nargin<12
    if printouttext, disp(' --- reference dispersion = 0 V, rerr disp H'); end;
    refdispersion=zeros(size(indBPM),2);
    [l,~,~]=atlinopt(rerr,0,indBPM);
    refdispersion(1,:)=arrayfun(@(a)a.Dispersion(1),l);
end


if scalefactor<0 || scalefactor>1
    if printouttext
        disp(' --- scale factor out of range. Set to 1.0'); end;
    scalefactor=1.0;
end


if correctflags(1) % dpp correction
    rmsel=[1 2 3 7 8 9];
else
    rmsel=[1 2 7 8];
end

% load or compute response matrix
if isempty(ModelRM)
    % get orbit RM
    if printouttext
        disp('get RM'); end;
    
    ModelRM=getresponsematrices(...
        rerr,...          %1 AT lattice
        indBPM,...      %2 bpm indexes in at lattice
        indHCor,...     %3 h cor indexes
        indVCor,...     %4 v cor indexes
        [],...     %5 skew cor indexes
        [],...     %6 quad cor indexes
        [],...
        inCOD,...       %7 initial coordinates
        rmsel...        %8 specifiy rm to be computed
        );
    
    if ~correctflags(1)
        
        ModelRM.OrbHDPP=[];
        ModelRM.OrbVDPP=[];
        ModelRM.DispHDPP=[];
        ModelRM.DispVDPP=[];
    end
    
end

% load RM computed by getresponsematrices

ormH=ModelRM.OrbHCor;
ormV=ModelRM.OrbVCor;
drmH=ModelRM.DispHCor;
drmV=ModelRM.DispVCor;
% kval=ModelRM.kval;
dppH=ModelRM.OrbHDPP;
dppV=ModelRM.OrbVDPP;
dppHd=ModelRM.DispHDPP;
dppVd=ModelRM.DispVDPP;
% delta=ModelRM.delta;
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
            

% get initial orbit
o=findorbit6Err(rerr,indBPM,inCOD);
ox0=o(1,:);
oy0=o(3,:);
d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
dx0=d(1,:);
dy0=d(3,:);

%rerr0=rerr;

% iterate correction
Niter=size(neigSteerer,1);
for iter=1:Niter
    
    if printouttext
        disp(['Disp. Free Steering iter ' num2str(iter,'%d, ') ...
            ' n-eig: ' num2str(neigSteerer(iter,:),'%d, ') ...
            ' alpha: ' num2str(wdisp,'%2.2f ')]);
    end
    
    % initial corrector strengths
    corh0=atgetfieldvalues(rerr,indHCor,'PolynomB',{1,1});
    corv0=atgetfieldvalues(rerr,indVCor,'PolynomA',{1,1});
    
    % get current orbit
    o=findorbit6Err(rerr,indBPM,inCOD);
    ox=o(1,:);
    oy=o(3,:);
    d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
    dx=d(1,:);
    dy=d(3,:);
    
    % subtract reference orbit
    ox=ox-reforbit(1,:);
    oy=oy-reforbit(2,:);
    % subtract reference dispersion
    dx=dx-refdispersion(1,:);
    dy=dy-refdispersion(2,:);
    
    % weigths between orbit and dispersion
    ox=ox*(1-wdisp);
    oy=oy*(1-wdisp);
    dx=dx*(wdisp);
    dy=dy*(wdisp);
    
    % build RMs
    if correctflags(1) && correctflags(2) % dpp and mean0
        RMH=[ [ormH{1}*(1-wdisp);drmH{1}*(wdisp);ones(size(indHCor))] [dppH'*(1-wdisp);dppHd'*(wdisp);0] ];
        RMV=[ [ormV{3}*(1-wdisp);drmV{3}*(wdisp);ones(size(indVCor))] [dppV'*(1-wdisp);dppVd'*(wdisp);0] ];
    elseif correctflags(1) && ~correctflags(2)% dpp no mean 0
        RMH=[ [ormH{1}*(1-wdisp);drmH{1}*(wdisp)] [dppH'*(1-wdisp);dppHd'*(wdisp)] ];
        RMV=[ [ormV{3}*(1-wdisp);drmV{3}*(wdisp)] [dppV'*(1-wdisp);dppVd'*(wdisp)] ];
    elseif ~correctflags(1) && correctflags(2) % mean0 no dpp
        RMH=[ormH{1}*(1-wdisp);drmH{1}*(wdisp);ones(size(indHCor))];
        RMV=[ormV{3}*(1-wdisp);drmV{3}*(wdisp);ones(size(indVCor))];
    elseif ~correctflags(1) && ~correctflags(2) % no dpp no mean0
        RMH=[ormH{1}*(1-wdisp);drmH{1}*(wdisp)];
        RMV=[ormV{3}*(1-wdisp);drmV{3}*(wdisp)];
    end
    
    % compute correction
    if correctflags(2) % mean 0
        dch=qemsvd_mod(RMH,[ox';dx';0],neigSteerer(1));
        dcv=qemsvd_mod(RMV,[oy';dy';0],neigSteerer(2));
    else % no constraint on correctors mean
        dch=qemsvd_mod(RMH,[ox';dx'],neigSteerer(1));
        dcv=qemsvd_mod(RMV,[oy';dy'],neigSteerer(2));
    end
    
    
    % get total correctors values and apply scaling
    if correctflags(1)
        hs=corh0-dch(1:end-1)*scalefactor;
        vs=corv0-dcv(1:end-1)*scalefactor;
        % energy deviation
        inCOD(5)=dch(end);
    else
        hs=corh0-dch*scalefactor;
        vs=corv0-dcv*scalefactor;
    end
    
    % limit steerers strengths
    if ~isempty(steererlimit)
        hs(abs(hs)>steererlimit(1))=steererlimit(1);
        vs(abs(vs)>steererlimit(2))=steererlimit(2);
    end
    
    % apply correction in lattice
    rcor=atsetfieldvalues(rerr,indHCor,'PolynomB',{1,1},hs);
    rcor=atsetfieldvalues(rcor,indVCor,'PolynomA',{1,1},vs);
    
    % lattice start point for next iteration
    rerr=rcor;
end

% get current orbit
o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);
d=finddispersion6Err(rcor,indBPM,indrfc,alpha,delta,inCOD);
dxc=d(1,:);
dyc=d(3,:);
    
Lh=atgetfieldvalues(rcor,indHCor,'Length');
Lv=atgetfieldvalues(rcor,indVCor,'Length');
hsL=hs.*Lh;
vsL=vs.*Lv;

if printouttext
    % display results
    disp(['      before' '    ' '-->' '    ' 'after'])
    disp(['oX: ' num2str(std(ox0-reforbit(1,:))*1e6,'%3.3f') ' -> ' num2str(std(oxc-reforbit(1,:))*1e6,'%3.3f') 'um']);
    disp(['oY: ' num2str(std(oy0-reforbit(2,:))*1e6,'%3.3f') ' -> ' num2str(std(oyc-reforbit(2,:))*1e6,'%3.3f') 'um']);
    disp(['dX: ' num2str(std(dx0-refdispersion(1,:))*1e3,'%3.3f') ' -> ' num2str(std(dxc-refdispersion(1,:))*1e3,'%3.3f') 'mm'])
    disp(['dY: ' num2str(std(dy0-refdispersion(2,:))*1e3,'%3.3f') ' -> ' num2str(std(dyc-refdispersion(2,:))*1e3,'%3.3f') 'mm'])
    disp(['    ' 'min' '    ' 'mean' '    ' 'max'])
    disp(['hs:'  num2str([min(hsL) mean(hsL) max(hsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['vs:'  num2str([min(vsL) mean(vsL) max(vsL)]*1e3,' %2.2f ') ' mrad'])
    disp(['dpp: ' num2str(inCOD(5))])
end


end
