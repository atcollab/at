function rerr=SetLargeErrorList(r0,seed,Nsig,factorerr,errnumber)
%
% function to set a given error list.
%
%see also: ApplyErrorsRand ApplyErrorsWave SetESRFAlgeAlignmentError

if nargin==4
    errnumber=1:8; %errors to apply from the list
end

if nargin==3
    factorerr=1; %multiply errors by this amount
    errnumber=[2:7,9]; %multiply errors by this amount
    
    disp('100% TDS list errors + current positions ALGE NO bpm err')
end

if factorerr==0
    factorerr=1e-12; % errors not sharp zero
end

if seed~=0
    disp(['Setting Random Stream to seed: ' num2str(seed)]);
    % set seed
    s = RandStream('mcg16807','Seed',seed);
    RandStream.setGlobalStream(s);
else
    disp('Using previously set random stream')
end

rerr=r0;

%% apply error wave

if find(errnumber==1)
    
    ie=1;
    
    wltouse=1:0.5:3;
    amplx=factorerr*0.6e-3;
    amplY=factorerr*0.6e-3;
    amplpsi=0*factorerr*0.6e-3;
    
    W=findspos(r0,length(r0)+1)./wltouse;
    
    A=amplx/length(W)*randn(size(W));
    errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
    errwavestruct(ie).type='x';
    errwavestruct(ie).A=A(end:-1:1);
    errwavestruct(ie).W=W;
    ie=ie+1;
    
    A=amplY/length(W)*randn(size(W));
    errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
    errwavestruct(ie).type='y';
    errwavestruct(ie).A=A(end:-1:1);
    errwavestruct(ie).W=W;
    ie=ie+1;
    
    A=amplpsi/length(W)*randn(size(W));
    errwavestruct(ie).indx=1:length(r0);%findcells(r0,'Class','Quadrupole');
    errwavestruct(ie).type='psi';
    errwavestruct(ie).A=A(end:-1:1);
    errwavestruct(ie).W=W;
    ie=ie+1;
    
    magindex=arrayfun(@(a)a.indx,errwavestruct,'un',0);
    type=arrayfun(@(a)a.type,errwavestruct,'un',0);
    A=arrayfun(@(a)a.A,errwavestruct,'un',0);
    W=arrayfun(@(a)a.W,errwavestruct,'un',0);
    
    rerr=ApplyErrorWave(...
        rerr,...
        magindex,...
        findcells(r0,'Class','Monitor'),...
        W,...
        A,...
        type);
    
end

%% define random errors structure
ie=1;
errstruct=[];

% %% GIRDERS (APPLY ALWAYS FIRST!)
% indg=findcells(r0,'FamName','GS');
% errstruct(ie).indx=indg;
% errstruct(ie).type='gx.gy';
% errstruct(ie).sigma=100*1e-6;
% ie=ie+1;
% errstruct(ie).indx=indg;
% errstruct(ie).type='gpsi';
% errstruct(ie).sigma=200*1e-6;
% ie=ie+1;

%% DIPOLES

if find(errnumber==2)
    
    % % DL
    indqm=find(atgetcells(r0,'FamName','DL\w*'));
    errstruct(ie).indx=indqm;
    errstruct(ie).type='x';
    errstruct(ie).sigma=100*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='y';
    errstruct(ie).sigma=100*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=200*1e-6;%500*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=indqm;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=1000*1e-6;
%     ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='dpb1';
    errstruct(ie).sigma=10*1e-4;
    ie=ie+1;
    
end

if find(errnumber==3)
    
    % DQ
    indqm=find(atgetcells(r0,'FamName','DQ\w*'));
    errstruct(ie).indx=indqm;
    errstruct(ie).type='x';
    errstruct(ie).sigma=70*1e-6;%70*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='y';
    errstruct(ie).sigma=50*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=150*1e-6;%200*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=indqm;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=500*1e-6;
%     ie=ie+1;
    
    errstruct(ie).indx=indqm;
    errstruct(ie).type='dpb1';
    errstruct(ie).sigma=10*1e-4;
    ie=ie+1;
    
    errstruct(ie).indx=indqm;
    errstruct(ie).type='dpb2';
    errstruct(ie).sigma=5*1e-4;
    ie=ie+1;
    
end

%% QUADRUPOLES

if find(errnumber==4)
    
    % moderate gradient quadrupoles
    indqm=find(atgetcells(r0,'FamName','Q[F-D][1-5]\w*'));
    errstruct(ie).indx=indqm;
    errstruct(ie).type='x';
    errstruct(ie).sigma=100*1e-6;%100*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='y';
    errstruct(ie).sigma=85*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=indqm;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=500*1e-6;
%     ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=150*1e-6;%200*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='dpb2';
    errstruct(ie).sigma=5*1e-4;
    ie=ie+1;
end

if find(errnumber==5)
    
    %high gradeint quadrupoles
    indqm=[find(atgetcells(r0,'FamName','Q[F-D][6-8]\w*'))];
    errstruct(ie).indx=indqm;
    errstruct(ie).type='x';
    errstruct(ie).sigma=50*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='y';
    errstruct(ie).sigma=70*1e-6;%70*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=indqm;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=500*1e-6;
%     ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=150*1e-6;%200*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indqm;
    errstruct(ie).type='dpb2';
    errstruct(ie).sigma=5*1e-4;
    ie=ie+1;
end

%% SEXTUPOLES

if find(errnumber==6)
    
    inds=find(atgetcells(r0,'Class','Sextupole'));
    errstruct(ie).indx=inds;
    errstruct(ie).type='x';
    errstruct(ie).sigma=70*1e-6;%70*1e-6;
    ie=ie+1;
    errstruct(ie).indx=inds;
    errstruct(ie).type='y';
    errstruct(ie).sigma=50*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=inds;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=1000*1e-6;
%     ie=ie+1;
    errstruct(ie).indx=inds;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=200*1e-6;%500*1e-6;
    ie=ie+1;
    errstruct(ie).indx=inds;
    errstruct(ie).type='dpb3';
    errstruct(ie).sigma=35*1e-4;
    ie=ie+1;
    
end

%% OCTUPOLES

if find(errnumber==7)
    
    indo=find(atgetcells(r0,'FamName','O[JF]\w*'));
    errstruct(ie).indx=indo;
    errstruct(ie).type='x.y';
    errstruct(ie).sigma=100*1e-6;
    ie=ie+1;
%     errstruct(ie).indx=indo;
%     errstruct(ie).type='s';
%     errstruct(ie).sigma=1000*1e-6;
%     ie=ie+1;
    errstruct(ie).indx=indo;
    errstruct(ie).type='psi';
    errstruct(ie).sigma=200*1e-6;%500*1e-6;
    ie=ie+1;
    errstruct(ie).indx=indo;
    errstruct(ie).type='dpb4';
    errstruct(ie).sigma=50*1e-4;
    ie=ie+1;
    
end

if find(errnumber==8)
    indm=find(atgetcells(r0,'Class','Monitor'));
    errstruct(ie).indx=indm;
    errstruct(ie).type='bpm';
    errstruct(ie).sigma=50*1e-6;
    ie=ie+1;
end

if find(errnumber==9)
  algedir='/mntdirect/_users/liuzzo/Matlab_Work/ATWORK/routines/esrfupgrade-gitrepo';
  algeactfile=fullfile(algedir,'Actual_Position_Simu.xlsx');
  rerr=SetESRFAlgeAlignmentError(rerr,algeactfile,'',seed);
end

if find(errnumber==10)
  algedir='/mntdirect/_users/liuzzo/Matlab_Work/ATWORK/routines/esrfupgrade-gitrepo';
  algeactfile=fullfile(algedir,'Nominal_Position_Simu.xlsx');
  rerr=SetESRFAlgeAlignmentError(rerr,algeactfile,'',seed);
end


if ~isempty(errstruct)
    
    %% set errors
    magindex=arrayfun(@(a)a.indx,errstruct,'un',0);
    type=arrayfun(@(a)a.type,errstruct,'un',0);
    sigma=arrayfun(@(a)a.sigma.*factorerr,errstruct,'un',0);
    
    rerr=ApplyErrorRand(...
        rerr,...
        magindex,...
        findcells(r0,'Class','Monitor'),...
        seed,...
        sigma,...
        Nsig,...
        type);
    
end




% figure;%('visible','off');
% atplot(rerr,@pltmisalignments);
% 
% labfig=['_seed' num2str(seed) '_Nsig' num2str(Nsig) '_scal' num2str(factorerr)];
% 
% saveas(gca,['Errors' labfig '.fig']);
% 
% try
%     export_fig(['Errors' labfig '.png'],'-transparent');
% catch
%     saveas(gca,['Errors' labfig '.png']);
% end

return