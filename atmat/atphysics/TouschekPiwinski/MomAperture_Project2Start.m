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
%       ETN: stability threshold for positive off energy particles
%       ETP: stability threshold for negative off energy particles
% Inputs:
%       THERING: ring used for tracking. Default: global THERING
%       REFPTS: REFPTS where to calculate the momentum acceptance.
%               Default 1:numel(THERING)+1;
%       nturns: Number of turns to track. Default 500
%       detole: resolution in energy acceptance. Default 1e-4
%       eu_ini: upper limit of the stability threshold. Default []
%               If not given it uses the linear energy acceptance from
%               ringpara
%       initcoord: [x y] starting transverse offset for the tracking.
%               Default [1e-6 1e-6]
%       verbose: boolean indicating verbose mode. Default false.
%       epsilon6D: if not passed, all particles are tracked.
%               If epsilon6D is given, we track for many turns only
%               the particles with 6D coordinates different by epsilon6D
%               after the 1st turn.
%
% Other functions in the file:
%
% Loste = Multiorigin_ringpass_islost
% Returns a boolean array: tells whether the particle launched at
% the reference point refpts with positive and negative energy offset is
% lost or not.

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
    fprintf('Particles differing by less than %.3e are considered similar.\n', ...
             epsilon6D);
end

np = numel(REFPTS);
if verbose
    fprintf('Using %d reference points.\n',np);
    fprintf('Tracking %d turns.\n',nturns);
    fprintf('Iteration stops when the energy step is below %.3f.\n',detole);
    fprintf('Using init coords %.3f %.3f um as transverse offsets.\n',1e6*initcoord(1), ...
             1e6*initcoord(2));
end

% initial bipartition settings
res=ringpara(THERING);
EACC=res.delta_max;
es_ini=0; % lower limit of the stability threshold
et_ini=EACC/2; % starting guess of the stability threshold
if isempty(eu_ini)
    eu_ini=EACC;
    if verbose
        fprintf('Using the rf bucket height as unstable energy limit.\n');
    end
end
if verbose
    fprintf('Unstable energy limit set at start to %.3f%%.\n',100*eu_ini);
end

% bipatition method for multiple points
orbit=findorbit6(THERING,REFPTS);

% positive/negative branch
etp= et_ini*ones(np,1);
eup= eu_ini*ones(np,1);
esp= es_ini*ones(np,1);
etn=-et_ini*ones(np,1);
eun=-eu_ini*ones(np,1);
esn=-es_ini*ones(np,1);
de=1;
iteration=0;
while de>detole && iteration < 100
    iteration=iteration+1;
    if verbose
        seconds_initial=datetime('now');
        fprintf('Boundary search, iteration %d ...\n',iteration);
    end
    % L is true for particles lost on the track
    L=Multiorigin_ringpass_islost(  THERING, ...
                                    REFPTS, ...
                                    etp, ...
                                    etn, ...
                                    orbit, ...
                                    nturns, ...
                                    initcoord, ...
                                    epsilon6D, ...
                                    verbose ...
                                 );
    % split in positive and negative side of energy offsets
    Lp = L(1:2:2*np);
    Ln = L(2:2:2*np);
    % split in stable (es), unstable (eu) and test (et) energy
    esp(Lp==0) = etp(Lp==0);
    eup(Lp~=0) = etp(Lp~=0);
    etp = (esp+eup)/2;
    esn(Ln==0) = etn(Ln==0);
    eun(Ln~=0) = etn(Ln~=0);
    etn = (esn+eun)/2;
    % define new energy step
    dep = max(abs(esp-eup));
    den = max(abs(esn-eun));
    de = max([dep,den]);
    if verbose
        elapsed_time=seconds(time(between(seconds_initial,datetime('now'))));
        fprintf('%1.3f seconds. Energy resolution is %1.3e and stops at %1.3e\n',elapsed_time,de,detole);
    end
end
end


function Loste=Multiorigin_ringpass_islost( THERING, ...
                                            refpts, ...
                                            ep, ...
                                            en, ...
                                            orbit, ...
                                            nturns, ...
                                            initcoord, ...
                                            epsilon6D, ...
                                            verbose ...
                                          )
% Loste=Multiorigin_ringpass_islost( THERING, ...
%                                    refpts, ...
%                                    ep, ...
%                                    en, ...
%                                    orbit, ...
%                                    nturns, ...
%                                    initcoord, ...
%                                    epsilon6D, ...
%                                    verbose ...
%                                  )
% Returns a boolean array: tells whether the particle launched at
% the reference point refpts with positive and negative energy offset is
% lost or not.
%   Inputs:
%       -THERING: cell array Lattice used for the traking.
%       -refpts: array [npossx1] with the elements number where to start the traking.
%       -ep: array [npossx1] with the positive energy deviation at each refpts.
%       -en: array [npossx1] with the negative energy deviation at each refpts.
%       -orbit: array [npossx6] with the orbit at each refpts.  
%       -nturns: number of full turns to track.
%       -initcoord: deviation from the 6D closed orbit, same dor each refpts.
%       -epsilon6D: minimum 6D distance between particles
%       -verbose: print additional info
%   Output:
%       -Loste : bool array [npossx2], True for lost particles.
%                every pair Loste(2*k-1,2*k) for k = 1,...
%                corresponds to the positive and negative energy offset.

nposs=numel(refpts);
Loste=ones(1,2*nposs);
Rin=zeros(6,2*nposs);
Rout=zeros(6,2*nposs);
length_epsilon6D = length(epsilon6D);
if length_epsilon6D == 1
    tinyoffset=epsilon6D;
end

% first track the remaining portion of the ring
for ii=1:nposs
    Line=THERING(refpts(ii):end);
    Rin(:,2*ii-1) = orbit(:,ii) + [initcoord(1) 0 initcoord(2) 0 ep(ii) 0.0]';
    Rin(:,2*ii  ) = orbit(:,ii) + [initcoord(1) 0 initcoord(2) 0 en(ii) 0.0]';
    [Rout(:,2*ii-1:2*ii),Loste(2*ii-1:2*ii)] = ringpass( ...
                                                        Line, ...
                                                        Rin(:,2*ii-1:2*ii) ...
                                                       );
end
nalive1stturn = length(Loste)-sum(Loste);
Ralive1stturn = Rout(:,Loste==0);
    
% use particles that have survived to the ring end,
% filter them if necessary, and track them
trackonly_mask = logical(1:length(Ralive1stturn));
similarparticles_index = [];
if length_epsilon6D == 1
    % search for non numerically similar particles
    DiffR = squeeze(std(repmat(Ralive1stturn,[1 1 nalive1stturn]) ...
                            - repmat(reshape( ...
                                            Ralive1stturn, ...
                                            [6 1 nalive1stturn] ...
                                            ), ...
                                            [1 nalive1stturn 1] ...
                                    ) ...
                        ) ...
                    );
    allposs = (1:nalive1stturn)'*ones(1,nalive1stturn);
    similarposs = max(allposs.*(DiffR<tinyoffset));
     % get mask of non numerically similar (100 x eps) particles
    [trackonly_mask, ~, similarparticles_index]=unique(similarposs);
    % Loste1=losses_in_multiturn(similarparticles_index);
    if verbose
        fprintf('Speed up when discarding similar particles %.3f%%\n', ...
                100*(1-length(trackonly_mask)/length(Ralive1stturn)));
    end
    
end
% track
% note: first call is only to reuse the lattice
ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 1e-6 0]',1);
[~, losses_in_multiturn] =ringpass( THERING, ...
                                    Ralive1stturn(:,trackonly_mask), ...
                                    nturns, ...
                                    'reuse' ...
                                  );
if length_epsilon6D == 1
    % copy losses result for numerically similar particles
    Lostpaux=losses_in_multiturn(similarparticles_index);
else
    Lostpaux = losses_in_multiturn;
end
% Now group with the particles that did not pass the first turn
Loste(Loste==0)=Lostpaux;
end