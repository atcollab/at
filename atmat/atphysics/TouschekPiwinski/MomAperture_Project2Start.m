function [etn, etp]=MomAperture_Project2Start(THERING, varargin)
% MOMAPERTURE_PROJECT2START calculates the local momentum aperture
%
% MOMAPERTURE_PROJECT2START is a Bipartition search of the negative and 
% positive stability thesholds in the 5th dimension (relative energy).
%  -The 6D closed orbit is taken into account.
%  -Particles launched at different REFPTS along the ring are first projected
%  to the ring last element so that all particles can be tracked together.
%
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns,detole)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns,detole,eu_ini)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns,detole,eu_ini,initcoord)
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns,detole,eu_ini,initcoord,verbose)
%          
%
% Output:
%       etn: stability threshold for positive off energy particles
%       etp: stability threshold for negative off energy particles
% Inputs:
%       THERING: ring used for tracking. Default: global THERING
%       REFPTS: REFPTS where to calculate the momentum acceptance.Default 1:numel(THERING)+1;
%       nturns: Number of turns to track. Default 500
%       detole: resolution in energy acceptance. Default 1e-4
%       eu_ini: upper limit of the stability threshold. Default [] (uses the linear energy acceptance from ringpara)
%       initcoord: [x y] starting transverse offset for the tracking. Default [1e-6 1e-6]
%       verbose: boolean indicating verbose mode. Default false.
%
% Other functions in the file:
%
% Loste=Multiorigin_ringpass_islost(poss,e,nturns):
% Returns an array of zeros and ones: tells whether the particle launched at
% position poss(ii) with energy e(ii) is lost or not after nturns turns.
% default input variables

% 2024may30 Z.Marti at ALBA CELLS

REFPTS=1:numel(THERING)+1;
nturns=500;
detole=0.0001;
eu_ini=[]; 
initcoord=[1e-6 1e-6];
verbose=false;

if nargin>=1
    THERING=varargin{1};
    REFPTS=1:numel(THERING)+1;
end
if nargin>=2
    REFPTS=varargin{2};
end
if nargin>=3
    nturns=varargin{3};
end
if nargin>=4
    detole=varargin{4};
end
if nargin>=5
    eu_ini=varargin{5};
end
if nargin>=6
    initcoord=varargin{6};
end
if nargin>=7
    verbose=varargin{7};
end


np = numel(REFPTS);

% initial bipartition settings
res=ringpara(THERING);
EACC=res.delta_max;
es_ini=0; % lower limit of the stability threshold
et_ini=EACC/2; % starting guess of the stability threshold
if isempty(eu_ini)
    eu_ini=EACC;
end


%% bipatition method for multiple points
orbit=findorbit6(THERING,REFPTS);

% positive branch
et=et_ini*ones(np,1);
eu=eu_ini*ones(np,1);
es=es_ini*ones(np,1);
de=1;
iteration=0;
while de>detole && iteration < 100
    if verbose
        seconds_initial=datetime('now');
        fprintf('Positive branch search, iteration %d...',iteration);
    end
    L=Multiorigin_ringpass_islost(THERING,REFPTS,et,orbit,nturns,initcoord);
    es(L==0)=et(L==0);
    eu(L~=0)=et(L~=0);
    et=(es+eu)/2;
    de=max(abs(es-eu));
    if verbose
        elapsed_time=seconds(time(between(seconds_initial,datetime('now'))));
        fprintf('Elapsed time is %1.3f seconds. Energy resolution is %1.3e\n',elapsed_time,de);
    end
    iteration=iteration+1;
end
etp=et;%+orbit(5,:)';
etp(etp<1e-6)=1e-6;

% negative branch
et=-et_ini*ones(np,1);
eu=-eu_ini*ones(np,1);
es=-es_ini*ones(np,1);
de=1;
iteration=0;
while de>detole && iteration < 100
    if verbose
        seconds_initial=datetime('now');
        fprintf('Negative branch search, iteration %d...',iteration);
    end
    L=Multiorigin_ringpass_islost(THERING,REFPTS,et,orbit,nturns,initcoord);
    es(L==0)=et(L==0);
    eu(L~=0)=et(L~=0);
    et=(es+eu)/2;
    de=max(abs(es-eu));
    if verbose
        elapsed_time=seconds(time(between(seconds_initial,datetime('now'))));
        fprintf('Elapsed time is %1.3f seconds. Energy resolution is %1.3e\n',elapsed_time,de);
    end
    iteration=iteration+1;
end
etn=et;%+orbit(5,:)';
etn(etn>-1e-6)=-1e-6;
end


function Loste=Multiorigin_ringpass_islost(THERING,refpts,e,orbit,nturns,initcoord)
% Returns an boolean array: tells whether the particle nposs launched at
% position poss(ii) with energy e(ii) is lost or not.
%
%   Inputs:
%       -THERING: cell array Lattice used for the traking.
%       -refpts: array [npossx1] with the elements number where to start the traking.
%       -e: array [npossx1] with the relative energy deviation at each refpts.
%       -orbit: array [npossx6] with the orbit at each refpts.  
%       -nturns: number of full turns to track.
%       -initcoord: deviation from the 6D closed orbit, same dor each refpts.

nposs=numel(refpts);
Loste=ones(nposs,1);
Rin=zeros(6,nposs);
Rout=zeros(6,nposs);
small=0.000001;
tiny=100*eps;

% first track the remaining portion of the ring
for ii=1:nposs
    Line=THERING(refpts(ii):end);
    Rin(:,ii)=orbit(:,ii)+[initcoord(1) 0 initcoord(2) 0 e(ii) 0.0]';
    Rout(:,ii)=ringpass(Line,Rin(:,ii));
    Loste(ii)=sign(sum(isnan(Rout(:,ii))));
end
ngood1=nposs-sum(Loste);
Rin1=Rout(:,Loste==0);
%Rin1(6,:)=0;
    
%Just track particles that have survived to the ring end 
if ngood1==1
    [~, LOSS1] =ringpass(THERING,Rin1,nturns);
    Loste(Loste==0)=LOSS1; 
elseif ngood1>1
    % search for non numerically similar (100 x eps) particles
    DiffR=squeeze(std(repmat(Rin1,[1 1 ngood1])-repmat(reshape(Rin1,[6 1 ngood1]),[1 ngood1 1])));
    allposs=(1:ngood1)'*ones(1,ngood1);
    similarposs=max(allposs.*(DiffR<tiny));
    [non_rep_poss,~, rep_index]=unique(similarposs);
    
    %track non numerically similar (100 x eps) particles
    ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 small 0]',1);
    [~, rep_LOSS1] =ringpass(THERING,Rin1(:,non_rep_poss),nturns,'reuse');

    % copy result for numerically similar (100 x eps) particles
    Loste1=rep_LOSS1(rep_index);
    
    % Now group with the particles that did not pass the first turn
    Loste(Loste==0)=Loste1;
end

end

