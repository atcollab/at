function [rcor,inCOD,qs,ss]=atcorrectdispersion(...
    rerr,...
    indBPM,...
    indQCor,...
    indSCor,...
    inCOD,...
    neigSteerer,...
    correctflags,...
    scalefactor,...
    ModelRM,...
    refdispersion,...
    correctorslimit,...
    printouttext)
% function [...
%    rcor,...           1) corrected lattice
%    inCOD,...          2) initial COD (dpp is stored here)
%    qs,ss...           3) required normal and skew quad. strengths (total)
%    ]=atcorrectdispersion(...
%     rerr,...          1) AT lattice to correct
%     indBPM,...        2) Nbx1 bpm indexes
%     indHCor,...       3) Nhx1 hor. cor indexes
%     indVCor,...       4) Nvx1 ver. cor indexes
%     inCOD,...         5) 6x1 initial COD guess
%     neigSteerer,...   6) 2xNiter eigenvectors for correction H and V at
%                          each iteration (default: [Nh/2 Nv/2])
%     correctflags,...  7) correct [dpp mean0](default: [true true])
%     scalefactor,...   8) scale factor to correction (default: 1.0)
%     ModelRM,...       9) ModelRM.Orb(H/V)Cor = 4x1 cell of orbit response mat.
%                          ModelRM.Orb(H/V)DPP = 6x1 array of orbit
%                          response to dpp
%                          if [] compute RM (default: [])
%     refdispersion,...10) 2xNbpm reference dispersion to correct to 
%                           (default rerr dispersion)
%     correctorslimit  11) 2x1 limit of steerers abs(steerer)<steererlimit
%                           (default: [], no limits)
%     printouttext     12) if 1 or true, display rms orbit
%     )
%
% features impelemented:
% limit correctors strengths
% ddp correction
% sum of steerers = 0
% 6D dispersion with BPM errors
% initial coordinates
% correction to reference dispersions refdispersion
% use atsetfieldvalues, atgetcells
%
%
%see also: qemsvd_mod finddispersion6Err getresponsematrices



% response matrix kicks
kval=1e-5;
delta=1e-3;

alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
f0=rerr{indrfc(1)}.Frequency;

% default arguments
if nargin<12
    printouttext=true;
end
if nargin<11
    correctorslimit=[];
end

if nargin<4
    if printouttext
        disp('get BPM and Correctors indexes'); end;
    indBPM=finc(atgetcells(rerr,'Class','Monitor'));
    indQCor=finc(atgetcells(rerr,'Class','Quadrupole'));
    indSCor=finc(atgetcells(rerr,'iscorS','S'));
end

if nargin<5
    inCOD=[0 0 0 0 0 0]';
end

if nargin<6
    neigSteerer=[length(indQCor) length(indSCor)]/2;
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
    if printouttext, disp(' --- computing orbit Response matrix'); end;
    ModelRM=[];
end

if nargin<10
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
    rmsel=[9 10 11];
else
    rmsel=[1 2];
end

% load or compute response matrix
if isempty(ModelRM)
    % get orbit RM
    if printouttext
        disp('get orbit RM'); end;
    
    ModelRM=getresponsematrices(...
        rerr,...          %1 AT lattice
        indBPM,...      %2 bpm indexes in at lattice
        [],...     %3 h cor indexes
        [],...     %4 v cor indexes
        indSCor,...     %5 skew  cor indexes
        indQCor,...     %6 quad cor indexes
        [],...
        inCOD,...       %7 initial coordinates
        rmsel...        %8 specifiy rm to be computed
        );
    
    if ~correctflags(1)
        
        ModelRM.DispHDPP=[];
        ModelRM.DispVDPP=[];
    end
    
end

drmH=ModelRM.DispQCor;
drmV=ModelRM.DispSCor;
% kval=ModelRM.kval;
dppH=ModelRM.DispHDPP;
dppV=ModelRM.DispVDPP;
% delta=ModelRM.delta;
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));
     
% get initial dispersion
d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
dx0=d(1,:);
dy0=d(3,:);

qs0=atgetfieldvalues(rerr,indQCor,'PolynomB',{1,2});

% iterate correction
Niter=size(neigSteerer,1);
for iter=1:Niter
    
    if printouttext
        disp(['Dispersion correction iter ' num2str(iter,'%d, ') ...
            'n-eig: ' num2str(neigSteerer(iter,:),'%d, ')]);
    end
    
    % initial corrector strengths
    corq0=atgetfieldvalues(rerr,indQCor,'PolynomB',{1,2});
    cors0=atgetfieldvalues(rerr,indSCor,'PolynomA',{1,2});
    
    % get current orbit
    d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
    dx=d(1,:);
    dy=d(3,:);
    
    % subtract reference orbit
    dx=dx-refdispersion(1,:);
    dy=dy-refdispersion(2,:);
    
    % build RMs
    if correctflags(1) && correctflags(2) % dpp and mean0
        RMH=[ [drmH{1};ones(size(indQCor))] [dppH';0] ];
        RMV=[ [drmV{3};ones(size(indSCor))] [dppV';0] ];
    elseif correctflags(1) && ~correctflags(2)% dpp no mean 0
        RMH=[ drmH{1} dppH' ];
        RMV=[ drmV{3} dppV' ];
    elseif ~correctflags(1) && correctflags(2) % mean0 no dpp
        RMH=[drmH{1};ones(size(indQCor))];
        RMV=[drmV{3};ones(size(indSCor))];
    elseif ~correctflags(1) && ~correctflags(2) % no dpp no mean0
        RMH=drmH{1};
        RMV=drmV{3};
    end
    
    % compute correction
    if correctflags(2) % mean 0
        dch=qemsvd_mod(RMH,[dx';0],neigSteerer(1));
        dcv=qemsvd_mod(RMV,[dy';0],neigSteerer(2));
    else % no constraint on correctors mean
        dch=qemsvd_mod(RMH,dx',neigSteerer(1));
        dcv=qemsvd_mod(RMV,dy',neigSteerer(2));
    end
    
    
    % get total correctors values and apply scaling
    if correctflags(1)
        qs=corq0-dch(1:end-1)*scalefactor;
        ss=cors0-dcv(1:end-1)*scalefactor;
        % energy deviation
        dd=-dch(end);
    else
        qs=corq0-dch*scalefactor;
        ss=cors0-dcv*scalefactor;
    end
    
    % limit correctors strengths
    if ~isempty(correctorslimit)
        qs(abs(qs)>correctorslimit(1))=correctorslimit(1);
        ss(abs(ss)>correctorslimit(2))=correctorslimit(2);
    end
    
    % apply correction in lattice
    rcor=atsetfieldvalues(rerr,indQCor,'PolynomB',{1,2},qs);
    rcor=atsetfieldvalues(rcor,indSCor,'PolynomA',{1,2},ss);
    
    if correctflags(1)
       
        rcor=atsetfieldvalues(rcor,indrfc,'Frequency',f0-alpha*(dd)*f0);
        
        if printouttext
            disp(['Delta RF : ' num2str(-alpha*(dd)*f0) ' Hz']);
        end
    end
    
    % lattice start point for next iteration
    rerr=rcor;
end

% get current orbit
d=finddispersion6Err(rcor,indBPM,indrfc,alpha,delta,inCOD);
dxc=d(1,:);
dyc=d(3,:);

if printouttext
    % display results
    disp(['before' '    ' '-->' '    ' 'after'])
    disp(['dX: ' num2str(std(dx0-refdispersion(1,:))*1e3,'%3.3f') ' -> ' num2str(std(dxc-refdispersion(1,:))*1e3,'%3.3f') 'mm'])
    disp(['dY: ' num2str(std(dy0-refdispersion(2,:))*1e3,'%3.3f') ' -> ' num2str(std(dyc-refdispersion(2,:))*1e3,'%3.3f') 'mm'])
    disp(['    ' 'min' '    ' 'mean' '    ' 'max'])
    disp(['hs:'  num2str([min(qs-qs0) mean(qs-qs0) max(qs-qs0)],' %2.4f ') ' 1/m2'])
    disp(['vs:'  num2str([min(ss) mean(ss) max(ss)],' %2.4f ') ' 1/m2'])
    disp(['dpp: ' num2str(inCOD(5))])

end
