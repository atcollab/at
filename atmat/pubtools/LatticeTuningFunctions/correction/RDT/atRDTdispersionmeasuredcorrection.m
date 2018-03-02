function [rcor,inCOD,qs,ss]=atRDTdispersionmeasuredcorrection(...
    rerr,...
    rfit,...
    r0,...
    indBPM,...
    indQCor,...
    indSCor,...
    inCOD,...
    neigSteerer,...
    correctflags,...
    scalefactor,...
    wdisp,...
    ModelRM,...
    steererlimit,...
    printouttext)
%ATRDTDISPERSIONMEASUREDCORRECTION Makes correction of dispersion based on
%RDTS
% function [...
%    rcor,...            1) corrected lattice
%    inCOD,...           2) initial COD (dpp is stored here)
%    hs,vs...            3) required steerers strengths (total)
%    ]=atdispersionfreesteering(...
%     rerr,...           1) AT lattice to correct, dispersion will be taken
%                           from this lattice
%     rfit,...           2) AT lattice with fitted errors
%                           optics will be taken from this lattice
%     r0, ...            3) 2xNbpm reference rdt to correct to
%     indBPM,...         4) Nbx1 bpm indexes       (default: monitor)
%     indQCor,...        5) Nqx1 quad. cor indexes (default: quadrupole)
%     indSCor,...        6) Nsx1 skew. cor indexes (default: sextupole)
%     inCOD,...          7) 6x1 initial COD guess  (default: 6x1 zero)
%     neigSteerer,...    8) 2xNiter eigenvectors for correction H and V at
%                          each iteration (default: [Nh/2 Nv/2])
%     correctflags,...   9) correct [ mean0](default: [ true])
%     scalefactor,...   10) scale factor to correction (default: 1.0)
%     [wdisph wtunes wdispv],...
%                       11) dispersion and tune weight:
%                           dispersionH*wdisph and orbith*(1-wdisph-wtune)
%                           dispersionV*wdispv and orbith*(1-wdispv)                          
%                           (default: 0.7 0.1 0.7)
%     ModelRM,...       12) ModelRM.Disp(N/S)Quad = 6x1 cell of dispersion 
%                           response mat. if [] compute RM (default: [])
%                           (default 0*2xNb, or from r0 if reftune is r0)
%     steererlimit      13) 2x1 limit of steerers abs(steerer)<steererlimit
%                           (default: [], no limits)
%     printouttext      14) if 1 or true, display rms orbit
%     )
%
% features impelemented:
% - limit correctors strengths
% - ddp correction
% - sum of steerers = 0
% - 6D orbit with BPM errors
% - initial coordinate
% - correction to reference rdt tune dispersion from r0 lattice
% - retrival of normal and skew quadrupole components also from alignment
%   errors and rotations
% - use atsetfieldvalues, atgetcells
%
%
% http://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.14.034002
%
%see also: qemsvd_mod findorbit6Err getresponsematrices



% response matrix kicks
%kval=1e-5;
delta=1e-3;

% default arguments
if nargin<14
    printouttext=true;
end
if nargin<13
    steererlimit=[];
end

if nargin<6
    if printouttext
        disp('get BPM and Correctors indexes'); end;
    indBPM=finc(atgetcells(rfit,'Class','Monitor'));
    indQCor=finc(atgetcells(rfit,'Class','Quadrupole'));
    indSCor=finc(atgetcells(rfit,'Class','Sextupole'));
end

if nargin<7
    inCOD=[0 0 0 0 0 0]';
end

if nargin<8
    neigSteerer=[length(indQCor) length(indSCor)]/2;
end

if nargin<9
    correctflags=true;
end

if nargin<10
    if printouttext
        disp(' --- scale set to 1.0'); end;
    scalefactor=1.0;
end

if nargin<11
    if printouttext, disp(' --- wdisph=0.7 wtune=0.1 wdispv=0.7'); end;
    wdisp=[.7 .1 .7];
end

if nargin<12
    if printouttext, disp(' --- computing orbit Response matrix'); end;
    ModelRM=[];
end


if scalefactor<0 || scalefactor>1
    if printouttext
        disp(' --- scale factor out of range. Set to 1.0'); end;
    scalefactor=1.0;
end


% load or compute response matrix
if isempty(ModelRM)
    % get orbit RM
    if printouttext
        disp('get RM'); end;
    
    
        ModelRM=getresponsematrices(...
            rfit,...          %1 AT lattice
            indBPM,...      %2 bpm indexes in at lattice
            [],...     %3 h cor indexes
            [],...     %4 v cor indexes
            indSCor,...     %5 skew cor indexes
            indQCor,...     %6 quad cor indexes
            [],...
            inCOD,...       %7 initial coordinates
            [10 11 12]...        %8 specifiy rm to be computed
            );
        
    
end

% load RM computed by getresponsematrices

drmQ=ModelRM.DispQCor;
drmS=ModelRM.DispSCor;


tuneQ=[ModelRM.TuneQCor{1};ModelRM.TuneQCor{2}];

% quad RDT RM
[~,~,ind]=EquivalentGradientsFromAlignments6D(rfit,inCOD);
indAllQuad=ind;
indAllSkew=ind;

%indAllQuad=[indQCor indSCor];
%indAllSkew=[indQCor indSCor];

[respqx,respqz]=qemrdtresp_mod(rfit,indBPM,indAllQuad);    % RDT response matrix assumes K=1
QL=atgetfieldvalues(rfit,indAllQuad,'Length');          % quadrupole lengths
QL(QL==0)=1;% thin lens magnets

% convert response from KL to K as for dispersion response matrix
% this is needed to merge the RM with the dispersion RM in the final
% computation of correction.
lengthsmat=repmat(QL',length(indBPM),1);
respqx=respqx.*lengthsmat;
respqz=respqz.*lengthsmat;

[~,qkcor]=ismember(indQCor,indAllQuad);
rdtQ=[...
    real(respqx(:,qkcor));...
    imag(respqx(:,qkcor));...
    real(respqz(:,qkcor));...
    imag(respqz(:,qkcor))];


% skew RDT RM
[respsx,respsz]=semrdtresp_mod(rfit,indBPM,indAllSkew);    % RDT response matrix assumes K=1
SL=atgetfieldvalues(rfit,indAllSkew,'Length');          % quadrupole lengths
SL(SL==0)=1;% thin lens magnets
lengthsmat=repmat(SL',length(indBPM),1);
respsx=respsx.*lengthsmat;
respsz=respsz.*lengthsmat;

[~,skcor]=ismember(indSCor,indAllSkew);
rdtS=[...
    real(respsx(:,skcor));...
    imag(respsx(:,skcor));...
    real(respsz(:,skcor));...
    imag(respsz(:,skcor))];


inCOD=[0 0 0 0 0 0]';
[l,t,~]=atlinopt(r0,0,indBPM);
refdispersion=zeros(2,length(indBPM));
refdispersion(1,:)=arrayfun(@(a)a.Dispersion(1),l);
refdispersion(2,:)=arrayfun(@(a)a.Dispersion(3),l);
reftune=t;

[KQnoer,KSnoer,~]=EquivalentGradientsFromAlignments6D(r0,inCOD);
%KQnoer=atgetfieldvalues(r0,indAllQuad,'PolynomB',{1,2});
%KSnoer=atgetfieldvalues(r0,indAllSkew,'PolynomA',{1,2});

fx=respqx*KQnoer;
fz=respqz*KQnoer;
rdtvecq=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';

fx=respsx*KSnoer;
fz=respsz*KSnoer;
rdtvecs=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';

refrdt(1,:)=rdtvecq;
refrdt(2,:)=rdtvecs;



% get rdt vectors to correct
[KQi,KSi,~]=EquivalentGradientsFromAlignments6D(rfit,inCOD);
%KQ=atgetfieldvalues(rerr,indAllQuad,'PolynomB',{1,2});
%KS=atgetfieldvalues(rerr,indAllSkew,'PolynomA',{1,2});

fx=respqx*KQi;
fz=respqz*KQi;
rq0=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';

fx=respsx*KSi;
fz=respsz*KSi;
rs0=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';


alpha=mcf(rfit);
indrfc=find(atgetcells(rfit,'Frequency'));

% get initial dispersion

d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
dx0=d(1,:);
dy0=d(3,:);

% get initial tune
[~,t0,~]=atlinopt(rerr,0,1);

%rerr0=rerr;
 qs0=atgetfieldvalues(rfit,indQCor,'PolynomB',{1,2});
 ss0=atgetfieldvalues(rfit,indSCor,'PolynomA',{1,2});
 qse0=atgetfieldvalues(rerr,indQCor,'PolynomB',{1,2});
 sse0=atgetfieldvalues(rerr,indSCor,'PolynomA',{1,2});
      
% iterate correction
Niter=size(neigSteerer,1);
for iter=1:Niter
    
    if printouttext
        disp(['RDT Disp. Tune Steering iter ' num2str(iter,'%d, ') ...
            ' n-eig: ' num2str(neigSteerer(iter,:),'%d, ') ...
            ' alpha: ' num2str(wdisp,'%2.2f ')]);
    end
    
    % initial corrector strengths
    corq0=atgetfieldvalues(rfit,indQCor,'PolynomB',{1,2});
    cors0=atgetfieldvalues(rfit,indSCor,'PolynomA',{1,2});
    corqe0=atgetfieldvalues(rerr,indQCor,'PolynomB',{1,2});
    corse0=atgetfieldvalues(rerr,indSCor,'PolynomA',{1,2});
    
    
    % get current rdt vectors to correct
    [KQe,KSe,~]=EquivalentGradientsFromAlignments6D(rfit,inCOD);
    %KQ=atgetfieldvalues(rerr,indAllQuad,'PolynomB',{1,2});
    %KS=atgetfieldvalues(rerr,indAllSkew,'PolynomA',{1,2});
    
    fx=respqx*KQe;
    fz=respqz*KQe;
    rq=[...
        real(fx);...
        imag(fx);...
        real(fz);...
        imag(fz)]';
    
    fx=respsx*KSe;
    fz=respsz*KSe;
    rs=[...
        real(fx);...
        imag(fx);...
        real(fz);...
        imag(fz)]';
    
    % get current dispersion
    d=finddispersion6Err(rerr,indBPM,indrfc,alpha,delta,inCOD);
    dx=d(1,:);
    dy=d(3,:);
    % get current tune
    [~,t,~]=atlinopt(rerr,0,1);
    
    
    % subtract reference orbit
    rq=rq-refrdt(1,:);
    rs=rs-refrdt(2,:);
    % subtract reference dispersion
    dx=dx-refdispersion(1,:);
    dy=dy-refdispersion(2,:);
    % subtract reference tune
    t=t-reftune;
    
    % weigths between RDT, tune and dispersion
    rq=rq*(1-wdisp(1)-wdisp(2));
    rs=rs*(1-wdisp(3));
    dx=dx*(wdisp(1));
    dy=dy*(wdisp(3));
    t=t*(wdisp(2));
    
    % build RMs
    if  correctflags(1) % mean0 no dpp
        RMQ=[rdtQ*(1-wdisp(1)-wdisp(2));drmQ{1}*(wdisp(1));tuneQ*(wdisp(2));ones(size(indQCor))];
        %RMQ=[rdtQ*(1-wdisp(1));drmQ{1}*(wdisp(1));ones(size(indQCor))];
        RMS=[rdtS*(1-wdisp(3));drmS{3}*(wdisp(3));ones(size(indSCor))];
    elseif ~correctflags(1) % no dpp no mean0
        RMQ=[rdtQ*(1-wdisp(1)-wdisp(2));drmQ{1}*(wdisp(1));tuneQ*(wdisp(2))];
        %RMQ=[rdtQ*(1-wdisp(1));drmQ{1}*(wdisp(1))];
        RMS=[rdtS*(1-wdisp(3));drmS{3}*(wdisp(3))];
    end
    
    % compute correction
    if correctflags(1) % mean = 0
        vecq=[rq dx t 0]';
        %vecq=[rq dx 0]';
        vecs=[rs dy 0]';
    else % no constraint on correctors mean
        vecq=[rq dx t]';
        %vecq=[rq dx]';
        vecs=[rs dy]';
    end
    
    dcq=qemsvd_mod(RMQ,vecq,neigSteerer(iter,1));
    dcs=qemsvd_mod(RMS,vecs,neigSteerer(iter,2));
    
    % get total correctors values and apply scaling
    
    qs=corq0-dcq*scalefactor;
    ss=cors0-dcs*scalefactor;
    qse=corqe0-dcq*scalefactor;
    sse=corse0-dcs*scalefactor;
    
    % limit steerers strengths
    if ~isempty(steererlimit)
        qs(abs(qs)>steererlimit(1))=steererlimit(1);
        ss(abs(ss)>steererlimit(2))=steererlimit(2);
        qse(abs(qse)>steererlimit(1))=steererlimit(1);
        sse(abs(sse)>steererlimit(2))=steererlimit(2);
    end
    
    % apply correction in lattice fitted and errors (for dispersion)
    rfit=atsetfieldvalues(rfit,indQCor,'PolynomB',{1,2},qs);
    rfit=atsetfieldvalues(rfit,indSCor,'PolynomA',{1,2},ss);
    
    rerr=atsetfieldvalues(rerr,indQCor,'PolynomB',{1,2},qse);
    rerr=atsetfieldvalues(rerr,indSCor,'PolynomA',{1,2},sse);
   
    % lattice corrected
    rcor=rfit;
end


% get current rdt vectors to correct
[KQ,KS,~]=EquivalentGradientsFromAlignments6D(rcor,inCOD);
%KQ=atgetfieldvalues(rcor,indQCor,'PolynomB',{1,2});
%KS=atgetfieldvalues(rcor,indAllSkew,'PolynomA',{1,2});

fx=respqx*KQ;
fz=respqz*KQ;
rqc=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';

fx=respsx*KS;
fz=respsz*KS;
rsc=[...
    real(fx);...
    imag(fx);...
    real(fz);...
    imag(fz)]';

% get current dispersion
d=finddispersion6Err(rcor,indBPM,indrfc,alpha,delta,inCOD);
dxc=d(1,:);
dyc=d(3,:);
% get current tune
[~,tc,~]=atlinopt(rcor,0,1);


if printouttext
    % display results
    disp(['        before' '    ' '-->' '    ' 'after'])
    disp(['rq: ' num2str(std(rq0-refrdt(1,:))*1e3,'%3.3f') ' -> ' num2str(std(rqc-refrdt(1,:))*1e3,'%3.3f') '']);
    disp(['rs: ' num2str(std(rs0-refrdt(2,:))*1e3,'%3.3f') ' -> ' num2str(std(rsc-refrdt(2,:))*1e3,'%3.3f') '']);
    disp(['dX: ' num2str(std(dx0-refdispersion(1,:))*1e3,'%3.3f') ' -> ' num2str(std(dxc-refdispersion(1,:))*1e3,'%3.3f') 'mm'])
    disp(['dY: ' num2str(std(dy0-refdispersion(2,:))*1e3,'%3.3f') ' -> ' num2str(std(dyc-refdispersion(2,:))*1e3,'%3.3f') 'mm'])
    disp(['tX: ' num2str((t0(1)-reftune(1)),'%3.3f') ' -> ' num2str((tc(1)-reftune(1)),'%3.3f') ''])
    disp(['tY: ' num2str((t0(2)-reftune(2)),'%3.3f') ' -> ' num2str((tc(2)-reftune(2)),'%3.3f') ''])
    disp(['    ' 'min' '    ' 'mean' '    ' 'max'])
    disp(['qs:'  num2str([min(qs-qs0) mean(qs-qs0) max(qs-qs0)]*1e0,' %2.2f ') ' 1/m2'])
    disp(['ss:'  num2str([min(ss-ss0) mean(ss-ss0) max(ss-ss0)]*1e0,' %2.2f ') ' 1/m2'])
    disp(['dpp: ' num2str(inCOD(5))])
    
   
%     figure;
%     subplot(4,1,1);
%     plot(rq0-refrdt(1,:),'r'); hold on;
%     plot(rqc-refrdt(1,:),'b'); 
%     legend('before','after')
%     subplot(4,1,2);
%     plot(rs0-refrdt(2,:),'r'); hold on;
%     plot(rsc-refrdt(2,:),'b'); 
%     subplot(4,1,3);
%     plot(dx0-refdispersion(1,:),'r'); hold on;
%     plot(dxc-refdispersion(1,:),'b'); 
%     subplot(4,1,4);
%     plot(dy0-refdispersion(2,:),'r'); hold on;
%     plot(dyc-refdispersion(2,:),'b'); 
%     saveas(gca,['RDTdispCor' num2str(wdisp,'%2.2f_') '.fig']);
%     export_fig(['RDTdispCor' num2str(wdisp,'%2.2f_') '.jpg']);
% %     
%     
end


end
