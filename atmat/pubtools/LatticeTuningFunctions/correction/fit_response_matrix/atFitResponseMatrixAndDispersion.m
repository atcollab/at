function [...
    rfit,...        % fitted lattice
    Kqn,...         % quad gradient errors
    Kqs,...         % skew gradient errors
    Kdh,...         % bending angle errors
    Kdv,...         % dipole rotation errors
    indquad,...     % index of Kqn and Kqs
    inddip...       % index of Kdh Kdv
    ]=atFitResponseMatrixAndDispersion(...
    rerr,...        1) lattice with errors to model
    r0,...          2) lattice without errors
    varargin)
%
% rfit =ATFITRESPONSEMATRIXANDDISPERSION(rerr,r0); computes an error model
%
% rfit models rerr using normal and skew quadrupole and dipole fields.
% this function can be copied and modified to adapt to fit any other
% quantity. The 'Numeric' derivative mode is then mandatory.
%
% INPUTS:
% rerr,...        1) lattice with errors to model
% r0,...          2) lattice without errors
% rerr and r0 must be the same lattice. only allowed difference are errors.
% if 'MeasuredRM' is provided, rerr is ignored.
% 'FitElement' field is added to the r0 structure to mark fit locations. if not present all quadrupole and bend are used
%
% OPTIONAL INPUTS (...,'NAME',value,... any order):
% inCOD         : zeros(6,1); guess input orbit
% indBPM        : bpm indexes in r0 (same as in rerr)
% indHCor       : hor cor. indexes in r0 (same as in rerr)
% indVCor       : hor cor. indexes in r0 (same as in rerr) for RM computation
% neig          : eigenvectors for fit [100, 100, 100] quad, dipole, skew errors
% nsubsets      : rm is divided in n subsets, fitted independently and averaged
% modecalc,     : 'Analytic'(default) or 'Numeric' (if function is modified to fit other parameters)
% speclab,...   : added to files saved with response derivative data
% measuredRM... : structure to store measured resonse matrix and dispersion
%                 see testMeasuredRMFit.m for an example of use.
%    measuredRM.rmfunctnorm : @fun(r,ib,ih,iv,~,txt) function handle to
%       compute a diagonal response matrix from model, used to compute RM derivative
%    measuredRM.rmfunctskew : @fun(r,ib,ih,iv,~,txt) function handle to
%       compute a off-diagonal response matrix from model, used to compute RM derivative
%    measuredRM.rmfunctdh   :@fun(r,ib,~,~,~,~) function handle to
%       compute a horizontal dispersion from model, used to compute derivative
%    measuredRM.rmfunctdv   :@fun(r,ib,~,~,~,~) function handle to
%       compute a vertical dispersion from model, used to compute derivative
%    measuredRM.RMH         : measured RM NBPMxNCor
%    measuredRM.RMV         : measured RM NBPMxNCor
%    measuredRM.DH          : measured hor dispersion
%    measuredRM.DV          : measured ver dispersion
%
%
% OUTPUTS:
%     rfit,...        % fitted lattice
%     Kqn,...         % quad gradient errors
%     Kqs,...         % skew gradient errors
%     Kdh,...         % bending angle errors
%     Kdv,...         % dipole rotation errors
%     indquad,...     % index of Kqn and Kqs
%     inddip...       % index of Kdh Kdv
%
%
%see also:


p = inputParser;

addRequired(p,'rerr',@iscell);
addRequired(p,'r0',@iscell);
addOptional(p,'inCOD',zeros(6,1),@isnumeric);
addOptional(p,'indBPM',find(atgetcells(r0,'Class','Monitor'))',@isnumeric);
addOptional(p,'indHCor',find(atgetcells(r0,'Class','Sextupole'))',@isnumeric);
addOptional(p,'indVCor',find(atgetcells(r0,'Class','Sextupole'))',@isnumeric);
addOptional(p,'neig',[100 100 100],@isnumeric);
addOptional(p,'nsubsets',1,@isnumeric);
addOptional(p,'modecalc','Analytic',@ischar);
addOptional(p,'speclab','',@ischar);
addOptional(p,'measuredRM',[]);
addOptional(p,'forceRMDerivativeRecomputation',false);

parse(p,rerr,r0,varargin{:});

rerr = p.Results.rerr;
r0 = p.Results.r0;
inCOD = p.Results.inCOD;
indBPM = p.Results.indBPM;
indHCor = p.Results.indHCor;
indVCor = p.Results.indVCor;
neig = p.Results.neig;
nsubsets = p.Results.nsubsets;
modecalc = p.Results.modecalc;
speclab = p.Results.speclab;
measuredRM = p.Results.measuredRM;
forceRMderivcalc = p.Results.forceRMDerivativeRecomputation;

%modecalc = 'Analytic' ; % 'Analytic' 'OAR' 'Numeric'
speclab = [speclab modecalc];

wdisph = 1;
wdispv = 1;
wtune = 1;

delta=1e-4;
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));

alpha0=mcf(r0);
f0=r0{indrfc(1)}.HarmNumber*PhysConstant.speed_of_light_in_vacuum.value/findspos(r0,length(r0)+1);


% "measure" RM, dispersion
if isempty(measuredRM)
    rmfunctnorm =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'normdisp',[wdisph wtune],txt);
    rmfunctskew =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'skewdisp',wdispv,txt);
    rmfunctdh =@(r,ib,~,~,~,~)getdisph6D(r,ib,indrfc,alpha,delta,inCOD);
    rmfunctdv =@(r,ib,~,~,~,~)getdispv6D(r,ib,indrfc,alpha,delta,inCOD);
    
    RMH=rmfunctnorm(rerr,indBPM,indHCor,indVCor,[],'computing hor. cor. RM');
    RMV=rmfunctskew(rerr,indBPM,indHCor,indVCor,[],'computing ver. cor. RM');
    DH=rmfunctdh(rerr,indBPM,[],[],[],'');
    DV=rmfunctdv(rerr,indBPM,[],[],[],'');
else
    rmfunctnorm = measuredRM.rmfunctnorm;
    rmfunctskew = measuredRM.rmfunctskew;
    rmfunctdh = measuredRM.rmfunctdh;
    rmfunctdv = measuredRM.rmfunctdv;
    
    RMH = measuredRM.RMH;
    RMV = measuredRM.RMV;
    DH = measuredRM.DH;
    DV = measuredRM.DV;
end

%% build "Response of RM to errors" for fit using OAR cluster.
indQuadsDeriv=find(atgetcells(r0,'Class','Quadrupole','Sextupole') & atgetcells(r0,'FitElement') )';
indDip=find(atgetcells(r0,'Class','Bend') & atgetcells(r0,'FitElement') )';
if isempty(indQuadsDeriv)
    warning('Empty quadrupole fit locations. using all quadrupoles');
    indQuadsDeriv=find(atgetcells(r0,'Class','Quadrupole','Sextupole') )';
end
if isempty(indDip)
    warning('Empty dipole fit locations. using all dipoles');
    indDip=find(atgetcells(r0,'Class','Bend'))';
end

indquad=indQuadsDeriv;
inddip=indDip;


% quad gradient error skew
QuadKn.index=indQuadsDeriv;
QuadKn.perturbation={'PolynomB',[1,2]};
QuadKn.delta=0;
QuadKn.dk=1e-6;

% quad gradient error skew
QuadKs.index=indQuadsDeriv;
QuadKs.perturbation={'PolynomA',[1,2]};
QuadKs.delta=0;
QuadKs.dk=1e-6;

% dipole angle
DipAng.index=indDip;
DipAng.perturbation={'BendingAngle',[1,1]};
DipAng.delta=0;
DipAng.dk=1e-6;

% dipole rotations
DipRot.index=indDip;
DipRot.perturbation=@(rr,ind,val)atsettilt(rr,ind,val); % Dipole rotation psiPAds(ipads).delta=0;
DipRot.delta=0;
DipRot.dk=1e-5;

    function CalcDerivatives(r0,forceRMderivcalc)
        % if files do not exist already
        if ~exist([speclab 'Norm' '.mat'],'file') || ...
                ~exist([speclab 'Skew' '.mat'],'file') || ...
                ~exist([speclab 'DH' '.mat'],'file') || ...
                ~exist([speclab 'DV' '.mat'],'file') || forceRMderivcalc
            
            switch modecalc
                
                case 'Analytic' % analytic RM derivative
                    
                    % warning('This option is available only for dipole and quadrupole errors fit (specify better !!!)');
                    
                    % get anlytical RM response
                    dpp = 0;
                    [...
                        dX_dq, dY_dq, ...
                        dDx_dq, ~, ~, ...
                        dXY_ds, dYX_ds, ...
                        ~, dDy_ds, ...
                        ~...
                        ]=AnalyticResponseMatrixDerivative(...
                        r0',dpp,...
                        indBPM',... % bpm
                        indHCor,...  % correctors
                        indquad,... % quadrupoles
                        indDip(1)',... % no dipoles, later
                        indquad);
                    
                    %  select relevant quadrupoles
                    
                    qind=find(atgetcells(r0,'Class','Quadrupole'))'; % response matrix computed always at all quadrupoles
                    quadforresponse=find(ismember(qind,indquad)); % quadrupoles to use for fit amoung all
                    nq = length(indquad);
                    dval = ones(size(indquad));
                   
                    % warning('dispersion/quadK not ok yet. never used in SR before 2018');
                    for iq=1:nq
                        
                        ib = quadforresponse(iq); % use only selected quadrupoles
                        rha=dX_dq(:,:,ib);
                        rva=dY_dq(:,:,ib);
                        dxdqa=[dDx_dq(:,ib)/(alpha0*f0);0;0]/dval(iq);% analytic, no tune dispersion derivative very different from
                        drespnorm(:,iq)=[rha(:);rva(:);dxdqa];
                        rhsa=dXY_ds(:,:,ib);
                        rvsa=dYX_ds(:,:,ib);
                        dxdqsa=[dDy_ds(:,ib)/(alpha0*f0)]/dval(iq);% analytic, no tune dispersion derivative very different from
                        drespskew(:,iq)=-[rvsa(:);rhsa(:);dxdqsa];
                        
                    end
                    
                    % dispersion RM
                    nb=length(indBPM);
                    alldip=find(atgetcells(r0,'BendingAngle'));
                    Ndip=numel(alldip);
                    dipsel=sort([1:Ndip]);
                    bndidx=sort(alldip(dipsel));% only long part
                    bang=atgetfieldvalues(r0,alldip,'BendingAngle');
                    
                    % make sure order is as in qempanel vector.
                    [~,b]=ismember(indDip,bndidx);
                    
                    % dispersion response varying dipole angle and rotation
                    [N,S,~,~]=dDxyDthetaDQ(r0,indBPM,alldip,'magmodel','thick');
                    % negative sign, probably due to dipole sign convention
                    N_sort= N(:,dipsel(b));
                    % qempanel response in [m/Hz] instead of [m/%]
                    N_sort=N_sort/(-alpha0*f0);
                    % in qempanel the response is for scaling factors, not deltas
                    % % N_sort=N_sort.*repmat(bang(dipsel(b))',nb,1);
                    
                    
                    S_sort=S(:,dipsel(b));
                    % qempanel response in [m/Hz] instead of [m/%]
                    S_sort=S_sort/(-alpha0*f0);
                    % in qempanel the response is for scaling factors, not deltas
                    S_sort=S_sort.*repmat(-bang(dipsel(b))',nb,1); % tilt not vertical angle
                    
                    rx = N_sort;
                    rz = S_sort;
                    
                    
                    
                    % save files to use with GenericRMFit
                    rmfunctnorm =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'normdisp',[wdisph wtune],txt);
                    rmfunctskew =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'skewdisp',wdispv,txt);
                    rmfunctdh =@(r,ib,~,~,~,~)getdisph6D(r,ib,indrfc,alpha,delta,inCOD);
                    rmfunctdv =@(r,ib,~,~,~,~)getdispv6D(r,ib,indrfc,alpha,delta,inCOD);
                    
                    rmcalcfun = rmfunctnorm;
                    PertArray = QuadKn;
                    dresp = drespnorm;
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'Norm.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    rmcalcfun = rmfunctskew;
                    PertArray = QuadKs;
                    dresp = drespskew;
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'Skew.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    rmcalcfun = rmfunctdh;
                    PertArray = DipAng;
                    dresp = rx ;
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'DH.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    rmcalcfun = rmfunctdv;
                    PertArray = DipRot;
                    dresp = rz ;
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'DV.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                case 'Numeric' % numeric quad, dip derivative
                    
                    % save files to use with GenericRMFit
                    rmcalcfun = rmfunctnorm;
                    PertArray = QuadKn;
                    dresp = GenericRMDerivative(r0,...  lattice
                        indBPM,... RM bpm indexes
                        indHCor,... RM corrector indexes
                        indVCor,...
                        PertArray,... function for perturbation
                        rmcalcfun...     function to compute RM (default (getfullrespmatrixvector))
                        );
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    
                    save([speclab 'Norm.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    
                    rmcalcfun = rmfunctskew;
                    PertArray = QuadKs;
                    dresp = GenericRMDerivative(r0,...  lattice
                        indBPM,... RM bpm indexes
                        indHCor,... RM corrector indexes
                        indVCor,...
                        PertArray,... function for perturbation
                        rmcalcfun...     function to compute RM (default (getfullrespmatrixvector))
                        );
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'Skew.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    rmcalcfun = rmfunctdh;
                    PertArray = DipAng;
                    dresp = GenericRMDerivative(r0,...  lattice
                        indBPM,... RM bpm indexes
                        indHCor,... RM corrector indexes
                        indVCor,...
                        PertArray,... function for perturbation
                        rmcalcfun...     function to compute RM (default (getfullrespmatrixvector))
                        );
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'DH.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
                    rmcalcfun = rmfunctdv;
                    PertArray = DipRot;
                    dresp = GenericRMDerivative(r0,...  lattice
                        indBPM,... RM bpm indexes
                        indHCor,... RM corrector indexes
                        indVCor,...
                        PertArray,... function for perturbation
                        rmcalcfun...     function to compute RM (default (getfullrespmatrixvector))
                        );
                    LongPertArray=arrayfun(@(a)StretchPertArray(a),PertArray,'un',0);
                    PertArray=[LongPertArray{:}];
                    save([speclab 'DV.mat'],'dresp','PertArray','rmcalcfun','indBPM','indHCor','indVCor');
                    
            end
            
            
        end
        
    end

%% fit measured RM
% implement average of errors fitted over several (4x4, instead of 16) subsets of correctors.
% define masks to fit subsets
N=nsubsets;

maskh=zeros(length(indBPM),length(indHCor));
maskv=zeros(length(indBPM),length(indVCor));

fitsubsetsH=zeros(N,length(RMH));
fitsubsetsV=zeros(N,length(RMV));
for ii=1:N
    maskhii=maskh;
    maskvii=maskv;
    maskhii(:,ii:N:end)=maskhii(:,ii:N:end)+1;
    maskvii(:,ii:N:end)=maskvii(:,ii:N:end)+1;
    fitsubsetsH(ii,:)=[maskhii(:);maskvii(:);ones(length(indBPM),1);1;1];
    fitsubsetsV(ii,:)=[maskhii(:);maskvii(:);ones(length(indBPM),1)];
end



rfit=r0; % model starts from no errors.
Kqn=zeros(size(indQuadsDeriv))';
Kqs=zeros(size(indQuadsDeriv))';
Kdh=zeros(size(indDip))';
Kdv=zeros(size(indDip))';


CalcDerivatives(rfit,true);
[~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),1,fitsubsetsH,false,false);  % fit normal
Kqn=Kqn+DKqn;
CalcDerivatives(rfit,true);
[~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1); % fit normal
Kdh=Kdh+DKdh;
% CalcDerivatives(rfit,true);
% [~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),1,fitsubsetsH,false,false);  % fit normal
% Kqn=Kqn+DKqn;
% CalcDerivatives(rfit,true);
% [~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1,2); % fit normal
% Kdh=Kdh+DKdh;
CalcDerivatives(rfit,true);
[~,rfit,DKqs,~,~]=GenericRMFit(RMV,rfit,[speclab 'Skew.mat'],neig(3),1,fitsubsetsV,false,false); % fit skew
Kqs=Kqs+DKqs;
CalcDerivatives(rfit,true);
[~,rfit,DKdv,~,~]=GenericRMFit(DV,rfit,[speclab 'DV.mat'],neig(2),1); % fit skew
Kdv=Kdv+DKdv;
CalcDerivatives(rfit,true);
[~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),1,fitsubsetsH,false,false);  % fit normal
Kqn=Kqn+DKqn;
CalcDerivatives(rfit,true);
[~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1); % fit normal
Kdh=Kdh+DKdh;


%% compare set errors to fitted errors
comparefitmod=false;
if comparefitmod
    [resperr]=rmfunctnorm(rerr,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctnorm(rfit,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['RMerr: ' num2str(std(resperr))],....
        ['RMfit: ' num2str(std(respfit))]);
    title('response matrix deviation normal');
    
    [resperr]=rmfunctskew(rerr,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctskew(rfit,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['RMerr: ' num2str(std(resperr))],....
        ['RMfit: ' num2str(std(respfit))]);
    title('response matrix deviation skew');
    
    
    [resperr]=rmfunctdh(rerr,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctdh(rfit,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['DHerr: ' num2str(std(resperr))],....
        ['DHfit: ' num2str(std(respfit))]);
    title('response matrix deviation h dispersion');
    
    
    [resperr]=rmfunctdv(rerr,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctdv(rfit,indBPM,indHCor,indVCor,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['DVerr: ' num2str(std(resperr))],....
        ['DVfit: ' num2str(std(respfit))]);
    title('response matrix deviation v dispersion');
    
    
    figure; atplot(rerr);
    figure; atplot(rfit);
end

end