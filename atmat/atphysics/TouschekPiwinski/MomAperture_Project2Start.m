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
% [ETN, ETP]=MOMAPERTURE_PROJECT2START(THERING,REFPTS,nturns,detole,eu_ini,initcoord,verbose,epsilon6D)
%          
%
% Output:
%       etn: stability threshold for positive off energy particles
%       etp: stability threshold for negative off energy particles
% Inputs:
%       THERING: ring used for tracking. Default: global THERING
%       REFPTS: REFPTS where to calculate the momentum acceptance.
%               Default 1:numel(THERING)+1;
%       nturns: Number of turns to track. Default 500
%       detole: resolution in energy acceptance. Default 1e-4
%       eu_ini: upper limit of the stability threshold. Default []
%               (uses the linear energy acceptance from ringpara)
%       initcoord: [x y] starting transverse offset for the tracking. Default [1e-6 1e-6]
%       verbose: boolean indicating verbose mode. Default false.
%       epsilon6D: if not passed, all particles are tracked. if epsilon6D
%               is given, we track for many turns only the particles with
%               6D coordinates different by epsilon6D after the 1st turn.
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
epsilon6D = [];

if nargin>=2
    REFPTS=varargin{1};
end
if nargin>=3
    nturns=varargin{2};
end
if nargin>=4
    detole=varargin{3};
end
if nargin>=5
    eu_ini=varargin{4};
end
if nargin>=6
    initcoord=varargin{5};
end
if nargin>=7
    verbose=varargin{6};
end
if nargin>=8
    epsilon6D=varargin{7};
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

% positive/negative branch
etp=et_ini*ones(np,1);
eup=eu_ini*ones(np,1);
esp=es_ini*ones(np,1);
etn=et_ini*ones(np,1);
eun=eu_ini*ones(np,1);
esn=es_ini*ones(np,1);
de=1;
iteration=0;
while de>detole && iteration < 100
    if verbose
        seconds_initial=datetime('now');
        fprintf('Positive and negative boundary search, iteration %d...',iteration);
    end
    L=Multiorigin_ringpass_islost(THERING,REFPTS,etp,etn,orbit,nturns,initcoord,epsilon6D);
    esp(L==0)=etp(L==0);
    eup(L~=0)=etp(L~=0);
    etp=(esp+eup)/2;
    esn(L==0)=etn(L==0);
    eun(L~=0)=etn(L~=0);
    etn=(esn+eun)/2;
    dep=max(abs(esp-eup));
    den=max(abs(esn-eun));
    de = max([dep,den]);
    if verbose
        elapsed_time=seconds(time(between(seconds_initial,datetime('now'))));
        fprintf('Elapsed time is %1.3f seconds. Energy resolution is %1.3e\n',elapsed_time,de);
    end
    iteration=iteration+1;
end
end


function Loste=Multiorigin_ringpass_islost(THERING,refpts,ep,en,orbit,nturns,initcoord,epsilon6D)
% Returns an boolean array: tells whether the particle nposs launched at
% position poss(ii) with energy e(ii) is lost or not.
%
%   Inputs:
%       -THERING: cell array Lattice used for the traking.
%       -refpts: array [npossx1] with the elements number where to start the traking.
%       -ep: array [npossx1] with the positive energy deviation at each refpts.
%       -en: array [npossx1] with the negative energy deviation at each refpts.
%       -orbit: array [npossx6] with the orbit at each refpts.  
%       -nturns: number of full turns to track.
%       -initcoord: deviation from the 6D closed orbit, same dor each refpts.
%       -epsilon6D: minimum 6D distance between particles

nposs=numel(refpts);
Loste=ones(1,2*nposs);
Rin=zeros(6,2*nposs);
Rout=zeros(6,2*nposs);
if len(epsilon6D) == 1
    tinyoffset=epsilon6D;
end

% first track the remaining portion of the ring
for ii=1:nposs
    Line=THERING(refpts(ii):end);
    Rin(:,ii)=orbit(:,ii)+[initcoord(1) 0 initcoord(2) 0 ep(ii) 0.0]';
    Rin(:,ii+1)=orbit(:,ii)+[initcoord(1) 0 initcoord(2) 0 en(ii) 0.0]';
    [Rout(:,ii:(ii+1)),Loste(ii:(ii+1))] = ringpass(Line,Rin(:,ii:(ii+1)));
end
nalive1stturn=2*nposs-sum(Loste);
Ralive1stturn=Rout(:,Loste==0);
    
% track particles that have survived to the ring end
% note: first call is only to reuse the lattice
ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 1e-6 0]',1);
if len(epsilon6D) == 1
    % search for non numerically similar (100 x eps) particles
    DiffR = squeeze(std(repmat(Ralive1stturn,[1 1 nalive1stturn]) ...
            - repmat(reshape(Ralive1stturn,[6 1 nalive1stturn]),[1 nalive1stturn 1])));
    allposs = (1:nalive1stturn)'*ones(1,nalive1stturn);
    similarposs=max(allposs.*(DiffR<tinyoffset));
    [non_rep_poss,~, rep_index]=unique(similarposs);
    
    %track non numerically similar (100 x eps) particles
    [~, rep_LOSS1] =ringpass(THERING,Ralive1stturn(:,non_rep_poss),nturns,'reuse');

    % copy result for numerically similar (100 x eps) particles
    Loste1=rep_LOSS1(rep_index);
    
    % Now group with the particles that did not pass the first turn
    Loste(Loste==0)=Loste1;
else
    %track non numerically similar (100 x eps) particles
    [~, rep_LOSS1] =ringpass(THERING,Ralive1stturn(:,non_rep_poss),nturns,'reuse');

end

end

