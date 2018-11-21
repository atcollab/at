function [etn etp]=MomAperture_Projecect2Start(varargin)
% function [etn etp]=MomAperture_Projecect2Start(THERING,positions,nvolt,detole,eu_ini,initcoord);
% Bipartition search of the negative and positive stability thesholds.
%  -The 6D closed orbit is taken into account.
%  -Particles launched at different positions along the ring are projectet
%  to the ring last element so that all particles can be tracked together.
%
% Output:
%         etn: stability threshold for positive off energy particles
%         etp: stability threshold for negative off energy particles
% Inputs:
%         positions: positions where to calculate the momentum acceptance
%         nvolt: Number of turns to track
%         detole: resolution in energy acceptance
%         eu_ini: upper limit of the stability threshold
%         initcoord: [x y] starting positions for the tracking
%
% Other functions in the file:
%
% Loste=Multiorigin_ringpass_islost(poss,e,nvolt):
% Returns an array of zeros and ones: tells wether the particle launched at
% position poss(ii) with energy e(ii) is lost or not after nvolt turns.

% default input variables
nvolt=500;
detole=0.0001;
eu_ini=0.1; 
initcoord=[1e-6 1e-6];

if nargin==0
    global THERING;
    positions=1:numel(THERING)+1;
elseif nargin==1
    THERING=varargin{1};
    nvolt=varargin{1};
    positions=1:numel(THERING)+1;
elseif nargin==2
    THERING=varargin{1};
    positions=varargin{2};
elseif nargin==3
    THERING=varargin{1};
    positions=varargin{2};
    nvolt=varargin{3};
elseif nargin==4
    THERING=varargin{1};
    positions=varargin{2};
    nvolt=varargin{3};
    detole=varargin{4};
elseif nargin==5
    THERING=varargin{1};
    positions=varargin{2};
    nvolt=varargin{3};
    detole=varargin{4};
    eu_ini=varargin{5};
elseif nargin==6
    THERING=varargin{1};
    positions=varargin{2};
    nvolt=varargin{3};
    detole=varargin{4};
    eu_ini=varargin{5};
    initcoord=varargin{6};
end

np = numel(positions);
ss= findspos(THERING, positions);


% initial bipartition settings
res=atsummary(THERING,'NoDisplay');
EACC=res.energyacceptance; % EACC is the linear energy acceptance calculated with atsummary
es_ini=0; % lower limit of the stability threshold
et_ini=EACC; % starting guess of the stability threshold


%% bipatition method for multiple points
[MRING, MS, orbit] = findm66(THERING,positions);

% positive branch
et=et_ini*ones(np,1);
eu=eu_ini*ones(np,1);
es=es_ini*ones(np,1);
de=1;
while de>detole
    L=Multiorigin_ringpass_islost(THERING,positions,et,orbit,nvolt,initcoord);
    es(L==0)=et(L==0);
    eu(L~=0)=et(L~=0);
    et=(es+eu)/2;
    de=max(abs(es-eu));
end
etp=et+orbit(5,:)';


% negative branch
et=-et_ini*ones(np,1);
eu=-eu_ini*ones(np,1);
es=-es_ini*ones(np,1);
de=1;
while de>detole
    L=Multiorigin_ringpass_islost(THERING,positions,et,orbit,nvolt,initcoord);
    es(L==0)=et(L==0);
    eu(L~=0)=et(L~=0);
    et=(es+eu)/2;
    de=max(abs(es-eu));
end
etn=et+orbit(5,:)';

end


function Loste=Multiorigin_ringpass_islost(THERING,poss,e,orbit,nvolt,initcoord)
% Returns an array of zeros and ones: tells wether the particle launched at
% position poss(ii) with energy e(ii) is lost or not.

nposs=numel(poss);
Loste=ones(nposs,1);
Rin=zeros(6,nposs);
small=0.000001;
tiny=0.0000000000000000001;

% frist track the remaining portion of the ring
for ii=1:nposs
    Line=THERING(poss(ii):end);
    Rin(:,ii)=ringpass(Line,orbit(:,ii)+[initcoord(1) 0 initcoord(2) 0 e(ii) 0.0]');
    Loste(ii)=sign(sum(isnan(Rin(:,ii))));
end
ngood1=nposs-sum(Loste);

%Just track particles that have surbived to the ring end 
if ngood1==1
    Rin1=Rin(:,Loste==0);
    [Rfin, LOSS1] =ringpass(THERING,Rin1,nvolt);
    Loste(Loste==0)=LOSS1;
    
elseif ngood1>1
    Rin1=Rin(:,Loste==0);

    % search for non similar particles
    DiffR=squeeze(std(repmat(Rin1,[1 1 ngood1])-repmat(reshape(Rin1,[6 1 ngood1]),[1 ngood1 1])));
    allposs=(1:ngood1)'*ones(1,ngood1);
    similarposs=max(allposs.*(DiffR<tiny));
    [non_rep_poss,~, rep_index]=unique(similarposs);
    
    %track non similar particles
    ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 small 0]',1);
    [Rfin, rep_LOSS1] =ringpass(THERING,Rin1(:,non_rep_poss),nvolt,'reuse');

    % copy result for similar particles
    Loste1=rep_LOSS1(rep_index);
    
    % Now group with the particles that did not pass the first turn
    Loste(Loste==0)=Loste1;
end

end

