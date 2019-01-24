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
    inCOD,...       3) guess for initial coordinates
    indBPM,...      4) bpm indexes for rm fit
    indHCor,...     5) h correctors for rm fit
    indVCor,...     6) v correctors for rm fit
    neig,...        7) # eigenvectors for fit [quad, dip, skew]
    nsubsets,...    8) # subsets for fit [quad, skew] errors=<fit(subsets)>
    modecalc,...    9) 'Analytic', 'Numeric', 'OAR'(ESRF)
    speclab,...     10) label to distinguish rm files and avoid recomputation if already existing
    measuredRM...   11) measured RM structure with function handles to simulate the same measurement
    )
%
% Function for RM fit. BASED ON ASD OAR CLUSTER @ ESRF
%
% FitElement is a field added to the rerr structure to mark fit locations.
% fit quadrupoles + sextupoles & FitElement
% fit dipoles & FitElement
%
%see also:

% 
% p=inputParser();
% 
% p = addRequired('r');


if nargin<9
    modecalc='Numeric';
end
if nargin<10
    speclab='';
end
if nargin<11
    measuredRM = [];
end

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
indquad=indQuadsDeriv;
inddip=indDip;

forceRMderivcalc = false;

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

% if files do not exist already
if ~exist([speclab 'Norm' '.mat'],'file') || ...
        ~exist([speclab 'Skew' '.mat'],'file') || ...
        ~exist([speclab 'DH' '.mat'],'file') || ...
        ~exist([speclab 'DV' '.mat'],'file') || forceRMderivcalc
    
    switch modecalc
        
        case 'OAR'
            
            oarspecstring='-l walltime=6 -q asd';
            
            tmpName=[tempname '.mat']; % 'latat.mat';%
            save(tmpName,'r0'); %r0 with initialized fields for errors and offsets.
            
            % run OAR cluster RM derivative
            tic;
            RMFoldNorm=RunRMGenericDerivArray(...
                tmpName,...
                'r0',...
                50,...
                [QuadKn ],...
                indHCor,...
                indVCor,...
                indBPM,...
                rmfunctnorm,... norm quad rm
                [speclab 'Norm'],...
                oarspecstring);            % dpp
            
            RMFoldSkew=RunRMGenericDerivArray(...
                tmpName,...
                'r0',...
                50,...
                [QuadKs ],...
                indHCor,...
                indVCor,...
                indBPM,...
                rmfunctskew,... % skew quad RM
                [speclab 'Skew'],...
                oarspecstring);            % dpp
            
            RMFoldDH=RunRMGenericDerivArray(...
                tmpName,...
                'r0',...
                10,...
                [ DipAng],...
                indHCor,...
                indVCor,...
                indBPM,...
                rmfunctdh,... norm quad rm
                [speclab 'DH'],...
                oarspecstring);            % dpp
            
            RMFoldDV=RunRMGenericDerivArray(...
                tmpName,...
                'r0',...
                10,...
                [ DipRot],...
                indHCor,...
                indVCor,...
                indBPM,...
                rmfunctdv,... % skew quad RM
                [speclab 'DV'],...
                oarspecstring);            % dpp
            
            
            disp('Waiting for all files to be appropriately saved');
            pause(10);
            % collect data
            CollectRMGenericDerivOARData(RMFoldNorm,...
                fullfile(RMFoldNorm,['../' speclab 'Norm.mat']));
            CollectRMGenericDerivOARData(RMFoldSkew,...
                fullfile(RMFoldSkew,['../' speclab 'Skew.mat']));
            
            CollectRMGenericDerivOARData(RMFoldDH,...
                fullfile(RMFoldDH,['../' speclab 'DH.mat']));
            CollectRMGenericDerivOARData(RMFoldDV,...
                fullfile(RMFoldDV,['../' speclab 'DV.mat']));
            % remove supporting folders and files
            
            rmdir(RMFoldNorm,'s');
            rmdir(RMFoldSkew,'s');
            rmdir(RMFoldDH,'s');
            rmdir(RMFoldDV,'s');
            delete(tmpName);
            toc;
            
            
        case 'Analytic' % analytic RM derivative
          
            warning('This option is available only for dipole and quadrupole errors fit (specify better !!!)');
            
            % get anlytical RM response
            dpp = 0;
            [...
                dX_dq, dY_dq, ...
                dDx_dq, dDx_db, Dx, ...
                dXY_ds, dYX_ds, ...
                dDx_ds, dDy_ds, ...
                dDy_da...
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
            
            warning('dispersion/quadK not ok yet. never used in SR before 2018');
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

            [N,S,~,~]=dDxyDtheta(r0,indBPM,alldip,'magmodel','thick');
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

% 
% if nargout==0
%     disp('only matrices for fit computed');
%     return
% end

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
%
% [~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),2,fitsubsetsH,true,false);  % fit normal
% Kqn=Kqn+DKqn;
% [~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1); % fit normal
% Kdh=Kdh+DKdh;
% [~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),2,fitsubsetsH,true,false);  % fit normal
% Kqn=Kqn+DKqn;
% [~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),2); % fit normal
% Kdh=Kdh+DKdh;
% [~,rfit,DKqs,~,~]=GenericRMFit(RMV,rfit,[speclab 'Skew.mat'],neig(3),2,fitsubsetsV,true,false); % fit skew
% Kqs=Kqs+DKqs;
% [~,rfit,DKdv,~,~]=GenericRMFit(DV,rfit,[speclab 'DV.mat'],neig(4),1); % fit skew
% Kdv=Kdv+DKdv;


[~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),1,fitsubsetsH,false,false);  % fit normal
Kqn=Kqn+DKqn;
[~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1); % fit normal
Kdh=Kdh+DKdh;
[~,rfit,DKqn,~,~]=GenericRMFit(RMH,rfit,[speclab 'Norm.mat'],neig(1),1,fitsubsetsH,false,false);  % fit normal
Kqn=Kqn+DKqn;
[~,rfit,DKdh,~,~]=GenericRMFit(DH,rfit,[speclab 'DH.mat'],neig(2),1); % fit normal
Kdh=Kdh+DKdh;

% % if ~exist([speclab 'FitQNorm' '.mat'],'file') || ...
% %         ~exist([speclab 'FitQSkew' '.mat'],'file') || ...
% %         ~exist([speclab 'FitQDH' '.mat'],'file') || ...
% %         ~exist([speclab 'FitQDV' '.mat'],'file')
%
%     oarspecstring='-l walltime=6 -q nice';
%
%     % quad gradient error skew
%     QuadKn.index=indQuadsDeriv;
%     QuadKn.perturbation={'PolynomB',[1,2]};
%     QuadKn.delta=0;
%     QuadKn.dk=1e-6;
%
%     % quad gradient error skew
%     QuadKs.index=indQuadsDeriv;
%     QuadKs.perturbation={'PolynomA',[1,2]};
%     QuadKs.delta=0;
%     QuadKs.dk=1e-6;
%
%     % dipole angle
%     DipAng.index=indDip;
%     DipAng.perturbation={'BendingAngle',[1,1]};
%     DipAng.delta=0;
%     DipAng.dk=1e-6;
%
%     % dipole rotations
%     DipRot.index=indDip;
%     DipRot.perturbation=@(rr,ind,val)atsettilt(rr,ind,val); % Dipole rotation psiPAds(ipads).delta=0;
%     DipRot.delta=0;
%     DipRot.dk=1e-5;
%
%     tmpName=[tempname '.mat']; % 'latat.mat';%
%     save(tmpName,'rfit'); %r0 with initialized fields for errors and offsets.
%
%     % run OAR cluster RM derivative
%     tic;
%     RMFoldNorm=RunRMGenericDerivArray(...
%         tmpName,...
%         'rfit',...
%         50,...
%         [QuadKn ],...
%         indhsel,...
%         indvsel,...
%         indbsel,...
%         rmfunctnorm,... norm quad rm
%         [speclab 'FitQNorm'],...
%         oarspecstring);            % dpp
%
%     RMFoldSkew=RunRMGenericDerivArray(...
%         tmpName,...
%         'rfit',...
%         50,...
%         [QuadKs ],...
%         indhsel,...
%         indvsel,...
%         indbsel,...
%         rmfunctskew,... % skew quad RM
%         [speclab 'FitQSkew'],...
%         oarspecstring);            % dpp
%
%     RMFoldDH=RunRMGenericDerivArray(...
%         tmpName,...
%         'rfit',...
%         10,...
%         [ DipAng],...
%         indhsel,...
%         indvsel,...
%         indbsel,...
%         rmfunctdh,... norm quad rm
%         [speclab 'FitQDH'],...
%         oarspecstring);            % dpp
%
%     RMFoldDV=RunRMGenericDerivArray(...
%         tmpName,...
%         'rfit',...
%         10,...
%         [ DipRot],...
%         indhsel,...
%         indvsel,...
%         indbsel,...
%         rmfunctdv,... % skew quad RM
%         [speclab 'FitQDV'],...
%         oarspecstring);            % dpp
%
%
%     disp('Waiting for all files to be appropriately saved');
%     pause(10);
%     % collect data
%     CollectRMGenericDerivOARData(RMFoldNorm,...
%         fullfile(RMFoldNorm,['../' speclab 'FitQNorm.mat']));
%     CollectRMGenericDerivOARData(RMFoldSkew,...
%         fullfile(RMFoldSkew,['../' speclab 'FitQSkew.mat']));
%
%     CollectRMGenericDerivOARData(RMFoldDH,...
%         fullfile(RMFoldDH,['../' speclab 'FitQDH.mat']));
%     CollectRMGenericDerivOARData(RMFoldDV,...
%         fullfile(RMFoldDV,['../' speclab 'FitQDV.mat']));
%     % remove supporting folders and files
%
%     rmdir(RMFoldNorm,'s');
%     rmdir(RMFoldSkew,'s');
%     rmdir(RMFoldDH,'s');
%     rmdir(RMFoldDV,'s');
%     delete(tmpName);
%     toc;
%
% %end

% [~,rfit,DKqs,~,~]=GenericRMFit(RMV,rfit,[speclab 'FitQSkew.mat'],neig(3),1,fitsubsetsV,false,false); % fit skew
% Kqs=Kqs+DKqs;
% [~,rfit,DKdv,~,~]=GenericRMFit(DV,rfit,[speclab 'FitQDV.mat'],neig(4),1); % fit skew
% Kdv=Kdv+DKdv;
% %[~,rfit,DKqs,~,~]=GenericRMFit(RMV,rfit,[speclab 'FitQSkew.mat'],neig(3),1,fitsubsetsV,false,false); % fit skew
% %Kqs=Kqs+DKqs;
[~,rfit,DKqs,~,~]=GenericRMFit(RMV,rfit,[speclab 'Skew.mat'],neig(3),1,fitsubsetsV,false,false); % fit skew
Kqs=Kqs+DKqs;
[~,rfit,DKdv,~,~]=GenericRMFit(DV,rfit,[speclab 'DV.mat'],neig(2),1); % fit skew
Kdv=Kdv+DKdv;


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