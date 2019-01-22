function [dresp,RM]=GenericRMDerivative(...
    r,...  lattice
    indBPM,... RM bpm indexes
    indHCor,... RM corrector indexes
    indVCor,...
    PertFunArray,... function for perturbation
    rmcalcfun...     function to compute RM (default (nothing))
    )
%function [dresp]=GenericRMDerivative(...
%     r,...  lattice
%     indBPM,... RM bpm indexes
%     indHCor,... RM corrector indexes
%     indVCor,...
%     PertFunArray,... function for perturbation
%     rmcalcfun...     function to compute RM (default (getfullrespmatrixvector))
%    )
%
% computes derivative of function rmcalcfun for a given set of 
% perturbations defined by PertFunArray. This derivative is used for the
% determination of an error model based on the measurement of rmcalcfun. 
% 
% Example rmcalcfun = function to compute ResponseMatrix, PertFunArray =
% quadrupole gradients 
% 
% INPUTS:
% r     : AT lattice
% indBPM: index of BPMs
% indHCor: index of Horizontal correctors
% indVCor: index of Vertical correctors
% PertFunArray:  
%   1) an array of structures as follows:
%       PertFunArray(1).index=ind1
%       PertFunArray(1).perturbation={'PolynomB',[1,2]}
%       PertFunArray(1).delta (default=0)
%       PertFunArray(1).dk    (default=1e-5)
%
%   2) a function handle :
%       PertFunArray(1).perturbation=@fun(r,ind,dk)
%       Example of function handle perturbations: 
%       1) PertFunArray(1).perturbation=...
%           @(r,ind,dk)setcellstruct(r,'PolynomB',...
%                     PertFunArray(1).index(ind),K0(ind)+dk,1,2)
%       2) PertFunArray(1).perturbation=@(~,~,~)disp('hello')
%   3) an array mixing 1) and 2).
%
% rmcalcfun : a function handle with signature 
%             rmcalcfun(r,indBPM,indHCor,indVCor,dpp,printouttext);
%
% OUTPUTS:
% dresp : (PertFunArray.perturbation(rmcalcfun) - rmcalcfun ) / PertFunArray.dk
% RM  (optional): structure with fields:    
%     RM.dresp=dresp;
%     RM.indHCor=indHCor;
%     RM.indVCor=indVCor;
%     RM.indBPM=indBPM;
%     RM.PertArray=PertFunArray;
%     RM.rmcalcfun=rmcalcfun;
% 
%see also: 

if nargin<6
    warning('redefine input parsing')
    rmcalcfun=@(r,~,~)disp('missing function for response : resp = rmcalcfun(r,indBPM,indHCor,indVCor,dpp,printouttext); '); % for future changes to this function, or input
end

% defaults
deltaDEF=0;
DKDEF=1e-5;

    % model response matrix on energy
    resp0=rmcalcfun(r,indBPM,indHCor,indVCor,deltaDEF,'model orm');
% try
% catch exc
%     disp('Something wrong in rmcalcfun!')
%     getReport(exc)
%     
%     rmcalcfun=@getfullrespmatrixvector; % for future changes to this function, or input
%     
%     resp0=rmcalcfun(r,indBPM,indHCor,indVCor,deltaDEF,'model orm');
% end

l=arrayfun(@(a)length(a.index),PertFunArray);

dresp=[];%zeros(length(resp0),sum(l));
NPert = length(PertFunArray);
% loop perturbations
for ipert=1:NPert
    
    PerturbFun=PertFunArray(ipert);
    
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
    else
        PerturbationFunction=PerturbFun.perturbation;
        K0=zeros(size(PerturbFun.index));
    end
    
    
    if ~isfield(PerturbFun,'dk')
        DK=DKDEF;
    else
        DK=PerturbFun.dk;
    end
    
    if ~isfield(PerturbFun,'delta')
        delta=deltaDEF;
    else
        delta=PerturbFun.delta;
    end
    
        
    NQ=length(PerturbFun.index);
    
    % initialize derivative
    drespPert=zeros(length(resp0),NQ);
    
    for iq=1:NQ
        tic;
   
        % change i-th parameter
        
        rdk=PerturbationFunction(r,PerturbFun.index(iq),K0(iq)+DK);
        
        if iscell(PerturbFun.perturbation)
            namepert=[PerturbFun.perturbation{1} '(' num2str(PerturbFun.perturbation{2},'%d,')  ')'];
        else
            namepert=func2str(PerturbFun.perturbation);
        end
        % recompute RM
        dispstringcount=[namepert ': ' num2str(iq,'%4d') '/' num2str(NQ,'%4d')];
        
        % recompute resp0, if delta is different in each curve
        resp0=rmcalcfun(r,indBPM,indHCor,indVCor,delta,'model orm');
        % compute rm with perturbation
        respdk=rmcalcfun(rdk,indBPM,indHCor,indVCor,delta,dispstringcount);
        
        %fprintf(repmat('\b',1,length(dispstringcount)+1));
        
        % compute derivative
        drespPert(:,iq)=(respdk-resp0)./DK;
        
        t1 = toc;
        disp(['expected time to wait: ' num2str(t1*(NQ-iq)/60) ' min']);
    end
    %PertFunRMDerivArray(ipert)=PerturbFun;
    %PertFunRMDerivArray(ipert).RMderiv=drespPert;
    
    dresp=[dresp,drespPert];

    
end

if nargout==2
    RM.dresp=dresp;
    RM.indHCor=indHCor;
    RM.indVCor=indVCor;
    RM.indBPM=indBPM;
    RM.PertArray=PertFunArray;
    RM.rmcalcfun=rmcalcfun;
end

return