% test errors and correction functions
close all
clear all
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions');

addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/correction/');
addpath('/mntdirect/_machfs/liuzzo/CODE/LatticeTuningFunctions/errors/');

load ESRFLattice.mat

%% get RM
speclab='OrbitESRF';


modelrmfile=fullfile(pwd,['RMmodel' speclab '.mat']);%

if ~exist([modelrmfile],'file')
    
    ModelRM...
        =getresponsematrices(...
        ring,...
        indBPM,...
        indHCor,...
        indVCor,...
        indSCor,...
        indQCor,...
        [],...
        [0 0 0 0 0 0]',...
        [1 2 3]);
    
    save([modelrmfile],'ModelRM');
else
    load([modelrmfile],'ModelRM');
end

% set errors
ind=find(atgetcells(ring,'Class','Quadrupole','Sextupole'));
dx=5e-6*randn(size(ind));
dy=5e-6*randn(size(ind));

rerr=atsetshift(ring,ind,dx,dy);

%% apply correction

% no steerers limit, no reference orbit
[rcor,inCOD,hs,vs]=atcorrectorbit(rerr,...
    indBPM,...
    indHCor,...
    indVCor,...
    [0 0 0 0 0 0]',...
    [50 50],...
    [false true],...
    1.0,...
    ModelRM,...
    zeros(2,length(indBPM)),...
    [],...
    true);

% steerers limited, increase eigenvectors number and changing RF frequency
[rcor,...           % corrected lattice
    inCOD,....      % initial orbit guess after correction
    hs,...          % total horizontal steerers strenghts
    vs....          % total vertical steerers strengths
    ]=atcorrectorbit(....
    rerr,...        % lattice to be corrected
    indBPM,...      % BPM indexes
    indHCor,...     % horizontal steerers indexes
    indVCor,...     % vertical steerers indexes
    inCOD,...       % input 6D closed orbit guess
    [...            % several correction iterations 
    [10 20];...     % with different number of eigenvectors 
    [30 40];...     % for horizontal and vertical plane
    [50 60];...     % <-- iter 3, use 50 eig hor., 60 eig ver.
    [70 70];...     % <-- iter 4, use 70 eig hor., 70 eig ver.
    [80 80];...     % <-- iter 5, use 80 eig hor., 80 eig ver.
    [97 96];...
    [97 96]...
    ],...
    [true true],... % [do dpp correction, keep average of correctors zero] 
    1.0,...         % scale factor for correction
    ModelRM,...     % response matrix, if [], compute it 
    zeros(2,length(indBPM)),... % reference orbit to correct to
    [0.5e-3 0.5e-3],... % sterrer strengths limits
    true);          % verbosity flag


o=findorbit6Err(rerr,indBPM,inCOD);
oxe=o(1,:);
oye=o(3,:);

o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);
oyc=o(3,:);

sBPM=findspos(rcor,indBPM);
figure;subplot(2,1,1);
plot(sBPM,oxe,'.-');hold on; plot(sBPM,oxc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('hor. COD');
subplot(2,1,2);
plot(sBPM,oye,'.-');hold on; plot(sBPM,oyc,'.-');
legend('before','after');
xlabel('s [m]');
ylabel('ver. COD');
saveas(gca,'OrbitCor.fig');
% export_fig('OrbitCor.jpg','-r300');


% plot output

figure;
subplot(2,1,1);bar(hs);ylabel('hor.')
subplot(2,1,2);bar(vs);ylabel('ver.')

inCOD(5)

rcor0=rcor;
