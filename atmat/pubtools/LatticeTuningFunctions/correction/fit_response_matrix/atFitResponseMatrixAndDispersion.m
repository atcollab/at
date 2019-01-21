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
    indbsel,...     3) bpm indexes for rm fit
    indhsel,...     4) h correctors for rm fit
    indvsel,...     5) v correctors for rm fit
    neig,...        6) # eigenvectors for fit [quad, dip, skew]
    nsubsets,...    7) # subsets for fit [quad, skew] errors=<fit(subsets)>
    speclab...      8) label to distinguish rm files and avoid recomputation if already existing
    )
%
% Function for RM fit. BASED ON ASD OAR CLUSTER @ ESRF
%
% FitElement is a field added to the rerr structure to mark fit locations.
% fit quadrupoles + sextupoles & FitElement
% fit dipoles & FitElement
%
%see also:

if nargin<7
    speclab='';
end

delta=1e-4;
alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));

% "measure" RM, dispersion
rmfunctnorm =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'normdisp',[1e-2 1e-2],txt);
rmfunctskew =@(r,ib,ih,iv,~,txt)simulaterespmatrixmeasurements(r,inCOD,ib,ih,iv,'skewdisp',1e-2,txt);
rmfunctdh =@(r,ib,~,~,~,~)getdisph6D(r,ib,indrfc,alpha,delta,inCOD);
rmfunctdv =@(r,ib,~,~,~,~)getdispv6D(r,ib,indrfc,alpha,delta,inCOD);

RMH=rmfunctnorm(rerr,indbsel,indhsel,indvsel,[],'computing hor. cor. RM');
RMV=rmfunctskew(rerr,indbsel,indhsel,indvsel,[],'computing ver. cor. RM');
DH=rmfunctdh(rerr,indbsel,[],[],[],'');
DV=rmfunctdv(rerr,indbsel,[],[],[],'');

%% build "Response of RM to errors" for fit using OAR cluster.
indQuadsDeriv=find(atgetcells(r0,'Class','Quadrupole','Sextupole') & atgetcells(r0,'FitElement') )';
indDip=find(atgetcells(r0,'Class','Bend') & atgetcells(r0,'FitElement') )';
indquad=indQuadsDeriv;
inddip=indDip;

AnalyticRMderiv = false;

forceRMderivcalc = false;

if ~AnalyticRMderiv
    if ~exist([speclab 'Norm' '.mat'],'file') || ...
            ~exist([speclab 'Skew' '.mat'],'file') || ...
            ~exist([speclab 'DH' '.mat'],'file') || ...
            ~exist([speclab 'DV' '.mat'],'file') || forceRMderivcalc
        
        oarspecstring='-l walltime=6 -q asd';
        
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
        
        tmpName=[tempname '.mat']; % 'latat.mat';%
        save(tmpName,'r0'); %r0 with initialized fields for errors and offsets.
        
        % run OAR cluster RM derivative
        tic;
        RMFoldNorm=RunRMGenericDerivArray(...
            tmpName,...
            'r0',...
            50,...
            [QuadKn ],...
            indhsel,...
            indvsel,...
            indbsel,...
            rmfunctnorm,... norm quad rm
            [speclab 'Norm'],...
            oarspecstring);            % dpp
        
        RMFoldSkew=RunRMGenericDerivArray(...
            tmpName,...
            'r0',...
            50,...
            [QuadKs ],...
            indhsel,...
            indvsel,...
            indbsel,...
            rmfunctskew,... % skew quad RM
            [speclab 'Skew'],...
            oarspecstring);            % dpp
        
        RMFoldDH=RunRMGenericDerivArray(...
            tmpName,...
            'r0',...
            10,...
            [ DipAng],...
            indhsel,...
            indvsel,...
            indbsel,...
            rmfunctdh,... norm quad rm
            [speclab 'DH'],...
            oarspecstring);            % dpp
        
        RMFoldDV=RunRMGenericDerivArray(...
            tmpName,...
            'r0',...
            10,...
            [ DipRot],...
            indhsel,...
            indvsel,...
            indbsel,...
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
        
    end
else % analytic RM derivative
    error('not working feature')
    % get anlytical RM response
[...
    dX_dq, dY_dq, ...
    dDx_dq, dDx_db, Dx, ...
    dXY_ds, dYX_ds, ...
    dDx_ds, dDy_ds...
    ]=CalcRespXXRespMat_thick_V2(...
    mach',dpp,...
    bpmidx',... % bpm
    hsidxlat,...  % correctors
    varidx,... % quadrupoles
    bndidx',...
    varidx);

% normalize 
[l,t,ch]=atlinopt(r,0,indHCor);
bx=arrayfun(@(a)a.beta(1),l);
by=arrayfun(@(a)a.beta(2),l);

o=findorbit6Err(r,indBPM,inCOD);
Ox=o(1,:);
Oy=o(3,:);
d=finddispersion6Err(r,indBPM,indrfc,alpha,delta,inCOD);
Dx=d(1,:);
Dy=d(3,:);

% get bpm resolution value 
bpmresx=atgetfieldvalues(r,indBPM,'Reading',{1,1});
bpmresy=atgetfieldvalues(r,indBPM,'Reading',{1,2});
LH=atgetfieldvalues(r,indHCor,'Length',{1,1});
LV=atgetfieldvalues(r,indVCor,'Length',{1,1});


OH=ormH{1}./repmat(kval./sqrt(bx).*LH',length(indBPM),1);%./repmat(bpmresx,length(indHCor),1)';
OV=ormV{3}./repmat(kval./sqrt(by).*LV',length(indBPM),1);%./repmat(bpmresy,length(indVCor),1)';
OHV=ormH{3}./repmat(kval./sqrt(bx).*LH',length(indBPM),1);%./repmat(bpmresy,length(indHCor),1)';
OVH=ormV{1}./repmat(kval./sqrt(by).*LV',length(indBPM),1);%./repmat(bpmresx,length(indVCor),1)';

respvectornorm=[ OH(:) ;... H orm
                 OV(:) ;... V orm
               ]; % response in a column vector

respvectorskew=[ OHV(:) ;... H orm
                 OVH(:) ;... V orm
                ]; % response in a column vector

rm=[OH, OVH;... H orm
    OHV, OV;... V orm
    ];
end

if nargout==0
    disp('only matrices for fit computed');
    return
end

%% fit measured RM
% implement average of errors fitted over several (4x4, instead of 16) subsets of correctors.
% define masks to fit subsets
N=nsubsets;

maskh=zeros(length(indbsel),length(indhsel));
maskv=zeros(length(indbsel),length(indvsel));

fitsubsetsH=zeros(N,length(RMH));
fitsubsetsV=zeros(N,length(RMV));
for ii=1:N
    maskhii=maskh;
    maskvii=maskv;
    maskhii(:,ii:N:end)=maskhii(:,ii:N:end)+1;
    maskvii(:,ii:N:end)=maskvii(:,ii:N:end)+1;
    fitsubsetsH(ii,:)=[maskhii(:);maskvii(:);ones(length(indbsel),1);1;1];
    fitsubsetsV(ii,:)=[maskhii(:);maskvii(:);ones(length(indbsel),1)];
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
    [resperr]=rmfunctnorm(rerr,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctnorm(rfit,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['RMerr: ' num2str(std(resperr))],....
        ['RMfit: ' num2str(std(respfit))]);
    title('response matrix deviation normal');
    
    [resperr]=rmfunctskew(rerr,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctskew(rfit,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['RMerr: ' num2str(std(resperr))],....
        ['RMfit: ' num2str(std(respfit))]);
    title('response matrix deviation skew');
    
    
    [resperr]=rmfunctdh(rerr,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctdh(rfit,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    
    figure;
    plot(resperr,'b');
    hold on;
    plot(respfit,'r');
    legend(['DHerr: ' num2str(std(resperr))],....
        ['DHfit: ' num2str(std(respfit))]);
    title('response matrix deviation h dispersion');
    
    
    [resperr]=rmfunctdv(rerr,indbsel,indhsel,indvsel,delta,...
        ['computed ORM with errors dpp: ' num2str(delta)]);
    [respfit]=rmfunctdv(rfit,indbsel,indhsel,indvsel,delta,...
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

return