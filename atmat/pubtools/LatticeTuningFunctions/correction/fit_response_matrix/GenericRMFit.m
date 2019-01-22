function [rcor,...  lattice with errors and fitted errors as correction
    rfit,...  lattice with fitted errors
    DK1,...   value of the gradient changes
    dr,...
    dr0...
    ]=GenericRMFit(...
    r,...  ( measured RM or lattice to simulate RM measurement
    r0,.... ( reference lattice, the model to fit, without any error)
    RM,... ( passed by load(filename) or second output of RMGenericDeriv)
    neigvectors,...  number of eigenvectors to use.
    totfititer,...   number of fit reiterations
    subsets,...  Nxlength(RM) masks for subsets
    erroraveragezero,...
    fixedchromaticity,...
    delta,...
    verbose...
    )
%
% fits errors in r0 to model r
%
% function [rcor,...  lattice with errors and fitted errors as correction
%     rfit,...  lattice with fitted errors
%     DK1,...   value of the gradient changes
%     dr,...
%     dr0...
%     ]=GenericRMFit(...
%     r,...             1) measured RM or lattice to simulate RM measurement
%     r0,....           2) reference lattice, the model to fit, without any error)
%     RM,...RMfile,...  3) variation of response matrix induced by errors
%                          RM.rmcalcfun: function handle to compute RM
%                          RM.PertArray: array of structures specifying the
%                                        lattice errors to be fitted. If 1
%                                        or more perturbation are
%                                        functions, totfititer is set to 1.
%                          RM.dresp    : RM variation with errors
%                          RM.indBPM   : monitors used in RM computation
%                          RM.indHCor  : h cor used in RM computation
%                          RM.indVCor  : v cor used in RM computation
%                          if passing filename RMfile, RM structure is saved in the file RMfile
%                          see RunRMGenericDerivArray, CollectRMGenericDerivOAR
%     neigvectors,...   4) number of eigenvectors to use for errors fit
%     totfititer,...    5) number of fit reiterations
%     subsets,...       6) Nxlength(RM) masks for subsets. fit parameters 
%                          are computed fitting each subsets and averaging.
%     erroraveragezero, 7) keep error average = 0
%     fixedchromaticity,8) keep fixed chromaticity (only if errors are
%                          sextupoles!)
%     delta...          9) delta p /p for RM off energy computation
%     verbose...        10) print out level 0, 1 (default), 2
%     )
%
% adds the errors defined in RM.PertArray to the AT lattice r0 to best fit
% the output of RM.rmcalcfun.
% the output rfit contains the specified errors and is as close as possible
% to the lattice r to be modelled.
%
%see also: CollectRMGenericDerivOAR RMGenericDeriv getfullrespmatrixvector

%erroraveragezero=false;
%fixedchromaticity=true;

tic;
if nargin<10
    verbose=1;
end
if nargin<9
    delta=0;
end
if nargin<8
erroraveragezero=false;
end
if nargin<7
fixedchromaticity=false;
end
if nargin<6
    subsets=[];
end

if nargin<5
    totfititer=5;
end

% switch for r being an AT lattice or a measured RM!
if isnumeric(r)
    respmeas=r;
    r=r0;
end

% load data from RM variable
if ischar(RM)
    RM=load(RM);
end

if ~isfield(RM,'rmcalcfun')
    rmfunct=@getfullrespmatrixvector;
else
    rmfunct=RM.rmcalcfun;
end

sel_indBPM=[];
sel_indHCor=[];
sel_indVCor=[];

dresp0=RM.dresp;
indBPM=RM.indBPM;
indHCor=RM.indHCor;
indVCor=RM.indVCor;

PertArray=RM.PertArray;

if isempty(subsets)
   subsets = ones(size(RM.dresp(:,1)))';
end

% if one or more perturbations are functions, set # fit iterations to 1
% this is due to the absolute setting of functional perturbations.
iscellperturb=arrayfun(@(a)iscell(a.perturbation),PertArray);

if find(~iscellperturb)
    totfititer=1;
end

% get fit point indexes
indCor=arrayfun(@(a)a.index,PertArray);
%[indFit,~,icFit]=unique(indCor,'stable');
indFit=indCor;

% display fitting param
if verbose>0
    cellpert=arrayfun(@(a)iscell(a.perturbation),PertArray);
    
    strtodispcell=arrayfun(@(a)(...
        ['fitting: ' r0{a.index(1)}.Class ': ' a.perturbation{1} '(' num2str(a.perturbation{2},'%d,') ')']...
        ),PertArray(cellpert),'un',0);
    
    strtodispfun=arrayfun(@(a)(...
        ['fitting: ' r0{a.index(1)}.Class ': ' func2str(a.perturbation)]...
        ),PertArray(~cellpert),'un',0);
    
    %disp(unique(strtodispcell)')
a=unique(strtodispcell)';
   
strfunshow=unique(strtodispfun)';
lstrfunsh=cellfun(@length,strfunshow);
strfunshow=cell2mat(cellfun(@(a)a(1:min(lstrfunsh)),strfunshow,'un',0));

%disp(strfunshow) 

try
    fitparammsg=[' ' a{:} ' ' strfunshow];
catch err
    fitparammsg='multiple parameters';
end

% fitparammsg={};
% for ia =1:length(a)
%     fitparammsg=[fitparammsg,[' ' a{ia} ' ' strfunshow(ia,:)]];
% end

end
% vector of indexes
%sk = cell2mat(arrayfun(@(x)find(indCor==x,1),indFit,'un',0));
%[~,sk]=ismember(indCor,indFit);

sk=1:length(indFit);

% fit sext and quads
rfit=r0;
rcor=r;

DK1=zeros(size(sk'));

%totfititer=10;
plotfittingproc=0;
%mode.neigs=100;

if isempty(sel_indBPM)
    sel_indBPM=indBPM;
end
if isempty(sel_indHCor)
    sel_indHCor=indHCor;
end
if isempty(sel_indVCor)
    sel_indVCor=indVCor;
end

% K10=getcellstruct(r0,'PolynomB',indCor,1,2); % K1 of quadruopoles in the fitted model
% K1err=getcellstruct(r,'PolynomB',indCor,1,2); % K1 of quadruopoles in the fitted model

% define elements of dresp to use in the correction.
[~,hind]=ismember(sel_indHCor,indHCor); % elements from RM to be used for correction
[~,vind]=ismember(sel_indVCor,indVCor); % elements from RM to be used for correction
[~,bind]=ismember(sel_indBPM,indBPM); % elements from RM to be used for correction


% compute orbit response matrix  MODEL
msg=[]; if verbose>1, msg=['computed ORM model dpp: ' num2str(delta)]; end
[resp0]=rmfunct(r0,sel_indBPM,sel_indHCor,sel_indVCor,delta,msg);

if ~exist('respmeas','var')
    % get orbit response matrix  WITH ERRORS
    msg=[]; if verbose>1, msg=['computed ORM  with errors dpp: ' num2str(delta)]; end
    [resperr]=rmfunct(r,sel_indBPM,sel_indHCor,sel_indVCor,delta,msg);
else
    if verbose>1, disp('Giving measured RM as input'); end
    resperr=respmeas;
    r=r0;
end


dr0=resperr-resp0;

% number of correctors and bpms
Nhc=length(sel_indHCor);
Nvc=length(sel_indVCor);
Nbpm=length(sel_indBPM);

if verbose>1
disp([' --- num eig for fit: ' num2str(neigvectors) '/' num2str(min(size(dresp0(:,sk)))) '  --- ']);
disp([' --- num param to fit: ' num2str(length(DK1)) '  --- ']);
end

%remove nan
indexnansresperr=isnan(resperr);

respfit=resp0;
% response matrix deviation to be corrected
dr=(resperr-respfit);

if fixedchromaticity
    [l,~,~]=atlinopt(r0,0,indFit);
    
    betax=arrayfun(@(a)a.beta(1),l);
    betay=arrayfun(@(a)a.beta(2),l);
    etax=arrayfun(@(a)a.Dispersion(1),l);
    etay=arrayfun(@(a)a.Dispersion(3),l);
end

for fitloopit=1:totfititer
    dr0=dr;
    
    % fit errors using SVD
    RM=dresp0(:,sk);
    vec=dr(:);
    if erroraveragezero
        RM=[RM; ones(size(sk)) ]; %
        vec=[vec; 0];
    end
    if fixedchromaticity
        RM=[RM; betax.*etax; betay.*etay]; %
        vec=[vec; 0; 0];
    end
    % fit subsets of response matrix and average results
    for isubset=1:size(subsets,1)
        
        sel=((subsets(isubset,:)==1)' & ~indexnansresperr );
        if erroraveragezero
            sel=[sel; true];
        end
        if fixedchromaticity
            sel=[sel; true; true];
        end
        
        % use backslah if all eigenvectors should be used (avoid SVD).
        if     neigvectors<min(size(RM(sel,:)))
            dkfit(:,isubset) = qemsvd_mod(RM(sel,:),vec(sel),neigvectors);
        else
            if verbose>1, disp('using \ operator'); end
            dkfit(:,isubset) = RM(sel,:) \ vec(sel) ;
        end
        
    end
    
    DK1=DK1+mean(dkfit,2);% absolute value for fit
   
   % [rfit,K1f,K1ff]=SetPerturbArray(r0,PertArray,+DK1(icFit));
   %  [rcor]=SetPerturbArray(r,PertArray,-DK1(icFit));
   
    [rfit,K1f,K1ff]=SetPerturbArray(r0,PertArray,+DK1);
    [rcor]=SetPerturbArray(r,PertArray,-DK1);
   
    
%     
%     if fitloopit~=1 && plotfittingproc
%         figure;
%         plot(K1err-K10,'k.-'); % set errors
%         hold on;
%         plot(K1f-K10,'rx-');% previous correction
%         plot(K1ff-K10,'bx-'); % current correction
%         plot(dkfit,'ro-');
%         
%         legend('set errors','previous iteration','fitted errors','fit values');
%         title('fitted gradient errors');
%         %saveas(gca,'SetFitErrorsAllerrAllQuadMat.fig');
%         %export_fig('SetFitErrorsAllerrAllQuadMat.png','-transparent');
%     end
    
    % get orbit response matrix  WITH Fitted Errors
    msg=[]; if verbose>1, msg=['computed ORM with fitted errors']; end
    [respfit]=rmfunct(rfit,sel_indBPM,sel_indHCor,sel_indVCor,delta,msg);
    dr=(resperr-respfit);
    
    if verbose>1,
        disp([fitparammsg ' : ' num2str(fitloopit) '/' num2str(totfititer)...
            ' residual rm: ' num2str(std(dr0(~indexnansresperr))) ...
            ' -> ' num2str(std(dr(~indexnansresperr))) ...
            ' rms cor: ' num2str(std(dkfit))]);
    elseif verbose>0,
        disp([fitparammsg ' : ' num2str(fitloopit) '/' num2str(totfititer)...
            ' residual rm: ' num2str(std(dr0(~indexnansresperr))) ...
            ' -> ' num2str(std(dr(~indexnansresperr))) ...
            ]);
    end
    
end


dr0=resperr(~indexnansresperr)-resp0(~indexnansresperr);
drfit=resperr(~indexnansresperr)-respfit(~indexnansresperr);


if  plotfittingproc
    [~,~,dxfit,dyfit]=getorbdispbeta(rfit,indBPM);
    [~,~,dxerr,dyerr]=getorbdispbeta(r,indBPM);
    [~,~,dx0,dy0]=getorbdispbeta(r0,indBPM);
    figure;
    sbpm=findspos(r,indBPM);
    plot(sbpm,dxerr-dx0,'ko--','DisplayName','errors');hold on;
    plot(sbpm,dxfit-dx0,'r+-','DisplayName','fitted');
    legend toggle;
    title('horizontal dispersion fitted vs errors')
    saveas(gca,'DispersionHFitvsErr.fig');
    export_fig('DispersionHFitvsErr.pdf','-transparent');
    
    figure;
    sbpm=findspos(r,indBPM);
    plot(sbpm,dyerr-dy0,'ko--','DisplayName','errors');hold on;
    plot(sbpm,dyfit-dy0,'r+-','DisplayName','fitted');
    legend toggle;
    title('vertical dispersion fitted vs errors')
    saveas(gca,'DispersionVFitvsErr.fig');
    export_fig('DispersionVFitvsErr.pdf','-transparent');
    
end
if verbose>1, toc; end
dr0=resperr-resp0;
drfit=resperr-respfit;

return


function [rdk,initialval,finalval]=SetPerturbArray(r,PertFunArray,DKfitval)
%
% apply perturbations set found from fit. vector as long as the sum of
% indexes in PertFunArray

rdk=r;
initialval=zeros(size(DKfitval));
finalval=zeros(size(DKfitval));
cursor=1;
% loop perturbations
for ipert=1:length(PertFunArray)
    
    PerturbFun=PertFunArray(ipert);
    
    dkpert=1;
    
    %     if isfield(PerturbFun,'dk')
    %         dkpert=PerturbFun.dk;
    %     else
    %         dkpert=1;
    %     end
    
    DKv=DKfitval(cursor:cursor+length(PerturbFun.index)-1).*dkpert;
    
    % switch pertrubFun may be a function or a cell array {fieldname value}
    if iscell(PerturbFun.perturbation)
        % model quadrupole gradients
        K0=getcellstruct(r,...
            PerturbFun.perturbation{1},...
            PerturbFun.index,...
            PerturbFun.perturbation{2}(1),...
            PerturbFun.perturbation{2}(2));
        
        
        PerturbationFunction=@(rr,ind,ddkk)setcellstruct(rr,...
            PerturbFun.perturbation{1},...
            ind,ddkk,...
            PerturbFun.perturbation{2}(1),...
            PerturbFun.perturbation{2}(2));
        
        rdk=PerturbationFunction(rdk,PerturbFun.index,K0+DKv);
        
    else
        PerturbationFunction=PerturbFun.perturbation;
        %K0=zeros(size(DKv));
        
        
        rdk=PerturbationFunction(rdk,PerturbFun.index,DKv);
    end
    
    
    initialval(cursor:cursor+length(PerturbFun.index)-1)=0;
    finalval(cursor:cursor+length(PerturbFun.index)-1)=DKv;
    cursor=cursor+length(PerturbFun.index);
end

return