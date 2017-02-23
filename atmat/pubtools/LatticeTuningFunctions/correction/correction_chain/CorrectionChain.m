function [...
    rcor,...            % corrected lattice
    ch,...              % final H cor values
    cv,...              % final V cor values
    cq,...              % final Quad cor values
    cs,...              % final Skew Quad cor values
    inCOD,...
    d0,de,dc...         % lattice data structures d0=no err, de=err, dc=cor
    ]=CorrectionChain(...
    rerr,...            %1  initial lattice
    r0,...              %2  model lattice
    indBPM,...          %3  bpm index
    indHCor,...         %4  h steerers index
    indVCor,...         %5  v steerers index
    indSkewQuadCor,...  %6  skew quad index
    indQuadCor,...      %7  quadrupole correctors index
    Neig,...            %8  number of eigen vectors [NeigorbitH, NeigorbitV, NeigQuadrdt, Neigdispv, Neigdisph,neig rdt corr, SkewQuadRDT]
    corrorder,...       %9  correction order 1: orbit, 2: tune, 3: skewquad disp v 4: quad disp h 5: quad RDT 6: skew RDT
    ModelRM,...          %10 response matrices
    speclab,...          %11 response matrices
    verbose)             %12 verbose (false): if true print out all relevat quantities after each step in corrorder
%
%
% performs a loop of corrections as described in corparam
%
% [...
%     rcor,...            % corrected lattice
%     ch,...              % final H cor values
%     cv,...              % final V cor values
%     cq,...              % final Quad cor values
%     cs...               % final Skew Quad cor valu
%     ]=CorrectionChain(...
%     rerr,...            %1  initial lattice
%     r0,...              %2  model lattice
%     indBPM,...          %3  bpm index
%     indHCor,...         %4  h steerers index
%     indVCor,...         %5  v steerers index
%     indSkewQuadCor,...  %6  skew quad index
%     indQuadCor,...      %7  quadrupole correctors index
%     Neig,...            %8  number of eigen vectors [NeigorbitH, NeigorbitV, NeigQuadrdt, Neigdispv, Neigdisph,neig rdt corr, SkewQuadRDT]
%     corrorder,...       %9  correction order 1: orbit, 2: tune, 3: skewquad disp v 4: quad disp h 5: quad RDT 6: skew RDT
%     ModelRM,...          %10 response matrices
%     speclab,...          %11 response matrices
%     verbose)             %12 verbose (false): if true print out all relevat quantities after each step in corrorder
% 
%
%
% number of eigenvalues for SVD
% neigSteererH     = Neig(1); % H orbit correction
% neigSteererV     = Neig(2); % V orbit correction
% neigSkew         = Neig(3); % skew quadrupoles vert. disp. correction
% neigQuad         = Neig(4); % quadrupoles hor. disp. correction, beta correction and phase correction
% neigQuadFit      = Neig(5); % quadrupole errors fit
% neigDipFit       = Neig(6); % dipole errors fit
% neigSkewFit      = Neig(7); % skew quad fit
%
%
% corrorder=[0:7];
%
%     '(-1 ): RF cavity frequency and time lag tuning '...
%     '( 0): open trajectory (finds closed orbit) '...
%     '( 1): orbit '...
%     '( 2): tune '...
%     '( 3): chromaticity '...
%     '( 4): dispersion '...
%     '( 5): dispersion free steering '...
%     '( 6): rdt + dispersion correction '...
%     '( 7): fit errors model and correct model quad RDT + dispersion (6) '
%
% if the correction fails, nothing is done to the lattice.
%
%see also:
% findrespmat
% qemsvd_mod
% atsetRFCavity
% atfirstturntrajectory
% atcorrectorbit
% fittunedelta2fam
% atmatchchromdelta
% atcorrectdispersion
% atdispersionfreesteering
% atRDTdispersioncorrection
% FitResponseMatrixAndDispersion
%


t0=tic;
disp('>>>>>  Correction sequence started  <<<<<')
Nbpm=length(indBPM);
NVcor=length(indVCor);
NHcor=length(indHCor);
NScor=length(indSkewQuadCor);
NQcor=length(indQuadCor);
disp(['     # BPM  : ' num2str(Nbpm) ]);
disp(['     # V cor: ' num2str(NVcor)]);
disp(['     # H cor: ' num2str(NHcor)]);
disp(['     # S cor: ' num2str(NScor)]);
disp(['     # Q cor: ' num2str(NQcor)]);
disp('>>>>>  --------------------------  <<<<<')

if nargin<12
    verbose=false;
end

if nargin<11
    speclab='';
end

if nargin<10
    ModelRM=[];
else
    disp(' --- - - - - - - - - - - - - - - - - - - ---')
    disp(' ---                                     ---')
    disp(' --- FAST MODE: MODEL RM ARE BEING USED! ---')
    disp(' ---                                     ---')
    disp(' --- - - - - - - - - - - - - - - - - - - ---')
end


if nargin<9
    
    corrorder=[0:3];
    % 0: open trajectory
    % 1: orbit
    % 2: tune
    % 3: chromaticity
    
end

if nargin<8  % default number of eigenvectors is 100
    disp('100 eigen vectors for all corrections')
    Neig = 100 * ones(7,1);
end
% number of eigenvalues for SVD
neigSteererH     =Neig(1); % H orbit correction
neigSteererV     =Neig(2); % V orbit correction
neigSkew         =Neig(3); % skew quadrupoles vert. disp. correction
neigQuad         =Neig(4); % quadrupoles hor. disp. correction, beta correction and phase correction
neigQuadFit      =Neig(5); % quadrupole errors fit
neigDipFit       =Neig(6); % dipole errors fit
neigSkewFit      =Neig(7); % skew quad fit

ch=atgetfieldvalues(rerr,indHCor,'PolynomB',{1,1});
cv=atgetfieldvalues(rerr,indVCor,'PolynomA',{1,1});
cq=atgetfieldvalues(rerr,indQuadCor,'PolynomB',{1,2});
cs=atgetfieldvalues(rerr,indSkewQuadCor,'PolynomA',{1,2});

rerrINIT=rerr; % initial lattice with errors to compute correction (PolynomB(2) stores everything!)

% display selected correction order
disp('>>>>>      Correction STRATEGY :       <<<<<')
disp('.')
for iorddisp=corrorder
    switch iorddisp
        case -1
            disp('          RF cavity')
        case 0
            disp('          open trajectory (steerers)')
        case 1
            disp('          orbit (steerers) ')
        case 2
            disp('          tune (quadrupoles, 2 families)')
        case 3
            disp('          chromaticity (sextupoles, 2 families)')
        case 4
            disp('          dispersion (quadrupoles)')
        case 5
            disp('          dispersion free steering (correctors)')
        case 6
            disp('          RDT + dispersion (quadrupoles) ')
        case 7
            disp('          Fit Quad+Dip Errors ')
            disp('          Correct RDT and Dispersion of fitted model ')
        otherwise
    end
end
disp('.')
disp('>>>>>  --------------------------  <<<<<')


inCOD=[0 0 0 0 0 0]';
[l,~,~]=atlinopt(r0,0,indBPM);
refdispersion=zeros(2,length(indBPM));
refdispersion(1,:)=arrayfun(@(a)a.Dispersion(1),l);
refdispersion(2,:)=arrayfun(@(a)a.Dispersion(3),l);

% perform correction
iicor=1;
for cor=corrorder
    tic;
    mesgcor=['Correction Step: ' num2str(iicor) '/' num2str(length(corrorder))];
    disp(mesgcor);
    iicor=iicor+1;
    
    rerr0=rerr;% for correction display.
    inCODe=inCOD;
    
    switch cor
        
        case -1
            %% set rfcavity
            
            % decide if radiation is on or off in the lattice by looking
            % for Rad in Pass Methods
            radon=any(cellfun(@(a)sum(ismember('Rad',a.PassMethod))==3,rerr));
            
            % get cavity settings
            indrfc=find(atgetcells(rerr,'Frequency'));
            rfv=sum(atgetfieldvalues(rerr,indrfc,'Voltage'));
            harm=atgetfieldvalues(rerr,indrfc(1),'HarmNumber');
            tlag=atgetfieldvalues(rerr,indrfc(1),'TimeLag');
            
            disp(['Set RF cavity. '....
                num2str(rfv*1e-6) ' MV, '...
                num2str(harm) ' buckets, '...
                num2str(radon) ' radiation '...
                ]);
            
            % set cavities rf frequency and time lag for lattice with
            % errors
          
            rerr=atsetRFCavityErr(rerr,rfv,radon,harm,inCOD);
            
            [...
                rerr,...
                inCOD...
                ]=atRFcorrection(...
                rerr,...
                indBPM,...
                inCOD,...
                [1 1 1 1 1],...
                1,...
                ModelRM);

        case 0
            %% OPEN TRAJECTORY
            excursion=1e-3;
            disp(['Open Trajectory Correction started. '....
                num2str(excursion*1e3) 'mm escursion accepted']);
            
            if isempty(ModelRM)
                ModelRMtr=r0;
            else
                ModelRMtr=ModelRM;
            end
            [rerr,inCOD]=atfirstturntrajectory(...
                rerr,...
                inCOD,...
                indBPM,...
                indHCor,...
                indVCor,...
                excursion,...
                200,...
                [false true],...
                ModelRMtr);
            
            ch=atgetfieldvalues(rerr,indHCor,'PolynomB',{1,1});
            cv=atgetfieldvalues(rerr,indVCor,'PolynomA',{1,1});
            
            try
            catch exc
                
                getReport(exc,'extended');
                error('Failed trajectory correction');
                
            end
            
        case 1
            %% ORBIT CORRECTION
            disp(['Steerers to fix orbit H ' ...
                ' using ' num2str(neigSteererH) ' eig']);
            disp(['Steerers to fix orbit V ' ...
                ' using ' num2str(neigSteererV) ' eig']);
            
            
            [rerr,inCOD]=atcorrectorbit(rerr,...
                indBPM,...
                indHCor,...
                indVCor,...
                inCOD,...
                [[floor(linspace(1,neigSteererH,10)),neigSteererH,neigSteererH];...
                [floor(linspace(1,neigSteererV,10)),neigSteererV,neigSteererV]]',...
                [false true],... dpp correction and mean to zero
                1.0,...
                ModelRM);
                        
            
            ch=atgetfieldvalues(rerr,indHCor,'PolynomB',{1,1});
            cv=atgetfieldvalues(rerr,indVCor,'PolynomA',{1,1});
            
            try   
            catch exc
                
                getReport(exc);
                warning('Failed Orbit correction')
                
            end
        case 2
            %% TUNE MATCH (2 families)
            
            try % tune rematch
                
                rerr=fittunedelta2fam(rerr,r0);
                
            catch exc
                
                getReport(exc);
                save('latticeFailedTunecor.mat','rerr','r0');
                disp('Could not match Tune');
                warning('Could not match Tune');
                
            end
            
        case 3
            %% chromaticity correction
            
            disp(' - - - -  chromaticty correction - - - - - ');
            disp('All SF and All SD moved by a constant value');
            
            try
                
                indS=find(atgetcells(r0,'Class','Sextupole'))';
                pbsxt=atgetfieldvalues(r0,indS,'PolynomB',{1,3});
                indSF=indS(pbsxt>0);
                indSD=indS(pbsxt<0);
                
                [~,~,chrom]=atlinopt(r0,0,1);disp(['Nominal chrom: ' num2str(chrom,'%2.3f, ')]);
                [~,~,chrome]=atlinopt(rerr,0,1);disp(['Initial chrom: ' num2str(chrome,'%2.3f, ')]);
                
                rerr=atmatchchromdelta(rerr,chrom,{indSF,indSD});
                
                [~,~,chromcor]=atlinopt(rerr,0,1);disp(['Final chrom: ' num2str(chromcor,'%2.3f, ')]);
                
                
            catch exc
                
                getReport(exc);
                warning('Failed chromaticty correction')
                
            end
            
        case 4
            %% DISPERSION CORRECTION
            disp(['Normal quadrupoles to fix dispersion H ' ...
                ' using ' num2str(neigQuad) ' eig']);
            disp(['Skew quadrupoles to fix dispersion V ' ...
                ' using ' num2str(neigSkew) ' eig']);
            
            
            [rerr,inCOD]=atcorrectdispersion(rerr,...
                indBPM,...
                indQuadCor,...
                indSkewQuadCor,...
                inCOD,...
                [[floor(linspace(20,neigQuad,5)),neigQuad,neigQuad];...
                [floor(linspace(20,neigSkew,5)),neigSkew,neigSkew]]',...
                [true true],... dpp correction and mean to zero
                1.0,...
                ModelRM,...
                refdispersion,...
                [],...
                true);
            
            
            cq=atgetfieldvalues(rerr,indQuadCor,'PolynomB',{1,2});
            cs=atgetfieldvalues(rerr,indSkewQuadCor,'PolynomA',{1,2});
            
            try
                
            catch exc
                
                getReport(exc);
                warning('Failed Orbit correction')
                
            end
            
        case 5
            %% ORBIT+DISPERSION CORRECTION (dispersion free steering)
            disp(['Steerers to fix orbit and dispersion H ' ...
                ' using ' num2str(neigSteererH) ' eig']);
            disp(['Steerers to fix orbit and dispersion V ' ...
                ' using ' num2str(neigSteererV) ' eig']);
            
            
            [rerr,inCOD]=atdispersionfreesteering(...
                rerr,...
                indBPM,...
                indHCor,...
                indVCor,...
                inCOD,...
                [[floor(linspace(20,neigSteererH,10)),neigSteererH,neigSteererH];...
                [floor(linspace(20,neigSteererV,10)),neigSteererV,neigSteererV]]',...
                [true true],...
                1.0,...
                0.9,...
                ModelRM,...
                zeros(2,length(indBPM)),...
                refdispersion,...
                [],...
                true);
            
            ch=atgetfieldvalues(rerr,indHCor,'PolynomB',{1,1});
            cv=atgetfieldvalues(rerr,indVCor,'PolynomA',{1,1});
            
            try
                
            catch exc
                
                getReport(exc);
                warning('Failed Orbit-dispersion (DFS) correction')
                
            end
        case 6
            %% RDT+DISPERSION CORRECTION
            disp(['Quadrupoles to fix RDT, tune and dispersion H ' ...
                ' using ' num2str(neigQuad) ' eig']);
            disp(['Steerers to fix RDT and dispersion V ' ...
                ' using ' num2str(neigSkew) ' eig']);
            
            [rerr,inCOD]=atRDTdispersioncorrection(...
                rerr,...
                r0,...
                indBPM,...
                indQuadCor,...
                indSkewQuadCor,...
                inCOD,...
                [[floor(linspace(1,neigQuad,5)),neigQuad,neigQuad];...
                [floor(linspace(1,neigSkew,5)),neigSkew,neigSkew]]',...
                [false],...
                1.0,...
                [0.8 0.1 0.8],...
                ModelRM);
            
            cq=atgetfieldvalues(rerr,indQuadCor,'PolynomB',{1,2});
            cs=atgetfieldvalues(rerr,indSkewQuadCor,'PolynomA',{1,2});
            
            try
                
            catch exc
                
                getReport(exc);
                warning('Failed RDT and dispersion correction')
                
            end
            
       
        case 7
            %% RDT+DISPERSION CORRECTION from lattice error model
                % fit lattice errors model
                [rfit]=FitResponseMatrixAndDispersion(...
                    rerr,...
                    r0,...
                    inCOD,...
                    indBPM,...
                    indHCor(1:9*2:end),... % 4 correctors, 1 every 8 cells
                    indHCor(1:9*2:end),...  % 4 correctors, 1 every 8 cells
                    [neigQuadFit,neigDipFit,neigSkewFit,neigDipFit],...
                    4,...
                    [speclab 'fitrm']);
                
                % get change of strength of correctors
                fq=atgetfieldvalues(rfit,indQuadCor,'PolynomB',{1,2});
                fs=atgetfieldvalues(rfit,indSkewQuadCor,'PolynomA',{1,2});
                
                % correct RDT and dispersion of fitted error model
                [~,inCOD,fcq,fcs]=atRDTdispersioncorrection(...
                    rfit,... <<--- fitted error model! not lattice with errors!
                    r0,...
                    indBPM,...
                    indQuadCor,...
                    indSkewQuadCor,...
                    inCOD,...
                    [[floor(linspace(1,neigQuad,5)),neigQuad,neigQuad];...
                    [floor(linspace(1,neigSkew,5)),neigSkew,neigSkew]]',...
                    [true],...
                    1.0,...
                    [0.8 0.1 0.8],...
                    ModelRM);
                
                %fcq=atgetfieldvalues(rfitcor,indQuadCor,'PolynomB',{1,2});
                %fcs=atgetfieldvalues(rfitcor,indSkewQuadCor,'PolynomA',{1,2});
                
                % store proposed correction
                dcq(1,:)=(fcq-fq);
                dcs(1,:)=(fcs-fs);
                
            
            % set delta correctors strength in lattice with errors.
            cq=atgetfieldvalues(rerr,indQuadCor,'PolynomB',{1,2});
            cs=atgetfieldvalues(rerr,indSkewQuadCor,'PolynomA',{1,2});
            
            cq=cq+dcq'; %  add proposed correction on fitted lattice
            cs=cs+dcs';
            rerr=atsetfieldvalues(rerr,indQuadCor,'PolynomB',{1,2},cq);
            rerr=atsetfieldvalues(rerr,indSkewQuadCor,'PolynomA',{1,2},cs);
            
            
            try     
            catch exc
                
                getReport(exc);
                warning('Failed model fit and RDT + dispersion correction')
                
            end
            
            
        
        otherwise
            warning([num2str(cor) ': not a possible correction. [0:7]: '...
                '( 0): open trajectory (finds closed orbit) '...
                '( 1): orbit '...
                '( 2): tune '...
                '( 3): chromaticity '...
                '( 4): dispersion '...
                '( 5): dispersion free steering '...
                '( 6): rdt + dispersion correction '...
                '( 7): fit errors model and correct model quad RDT + dispersion (6) '])
    end
    
    %     %% set corrector multipoles
    %     rerr=SetCorMult(rerr);
    %
    
    % %%apply correctors PS limits
    %rerr=SetCorLimits(rerr);
    if verbose
        %% display correction effect at every step
        DisplayCorrectionEffect(r0,rerr0,rerr,inCODe,inCOD,1:length(r0),indHCor,indVCor,indQuadCor,indSkewQuadCor);
    end
    
    disp(['Finished: ' mesgcor])
    toc;
end

d0=[];de=[];dc=[];
if ~verbose
    %% display correction effect from begining to end.
 
    [d0,de,dc]=DisplayCorrectionEffect(r0,rerrINIT,rerr,inCODe,inCOD,[1:length(r0)]',indHCor,indVCor,indQuadCor,indSkewQuadCor);
end
rcor=rerr;

tend=toc(t0);
disp(['Time for correction chain: ' num2str((tend-t0)/60) ' minutes'])
return