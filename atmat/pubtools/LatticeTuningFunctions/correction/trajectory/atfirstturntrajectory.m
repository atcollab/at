function [rclosed,inCOD]=atfirstturntrajectory(...
    ropen,...
    inCOD,...
    indBPM,...
    indHCor,...
    indVCor,...
    lim,...
    maxcornum,...
    correctflags,...
    ModelRM,...
    reforbit,...
    steererlimit,...
    printouttext)
% [
% rclosed,...           1) AT lattice with closed trajectory
% inCOD...              2) 6x1 input coordinate updated
% ]=atfirstturntrajectory(...
%     ropen,...         1) AT lattice
%     inCOD,...         2) 6x1 input coordinate
%     indBPM,...        3) bpm indexes
%     indHCor,...       4) h. cor indexed
%     indVCor,...       5) v. cor indexed
%     lim,...           6) maximum bpm reading (lim+2e-3m for extra search)
%                          (default: 4e-3m)
%     maxcornum,...     7) maximum correctors to use for correction
%                          (only last maxcornum correctors and BPMs are used)
%     correctflags,...  8) correct [dpp mean0](default: [true true])
%     ModelRM,...(r0)   9) output of getresponsematrice
%                          (no default. if no rm has been computed, this
%                          argument can be r0, an AT lattice cell array of
%                          structures without errors to use for
%                          computation of the RM.
%     reforbit,...     10) reference trajectory 2xnbpm (default:0)
%     steererlimit,... 11) 2x1 limit of steerers abs(steerer)<steererlimit
%                          (default: [], no limits)
%     printouttext     12) print text ( default: true)
%     )
%
% finds closed trajectory:
% correct trajectory with available BPM and correctors
% if not found, increase amplitude limit for bpm reading
% if still not found scan injection point for maximum number of turns
% if all bpm see signal, close orbit using the last 2 correctors in the
% lattice.
%
%see also: findrespmat findtrajectory6Err getresponsematrices

% features to implement
% average correctors to zero
% dpp correction
% steerers limits
% reference orbit
% limit number of bpm/correctors to use default to all

% features implemented
% actions if stuck at a given bpm
%       - scan for optimal inCOD
%       - increase limit of readable bpm signal
%       - reduce number of eigenvectors
%

% defaults
if nargin<12
    printouttext=true;
end

if nargin<2
    inCOD=[0 0 0 0 0 0]';
end
if nargin<3
    if printouttext
        disp('get BPM and Correctors indexes');end;
    indBPM=find(atgetcells(ropen,'Class','Monitor'));
end
if nargin<4
    indHCor=find(atgetcells(ropen,'iscorH','H'));
end
if nargin<5
    indVCor=find(atgetcells(ropen,'iscorV','V'));
end
if nargin<6
    lim=4e-3;
end
if nargin<7
    maxcornum=length(indBPM);
end
if nargin<8
    correctflags=[true true];
end

if nargin<10
    reforbit=zeros(2,length(indBPM));
end

if nargin<11
    steererlimit=[];
end

if iscell(ModelRM)
    r0=ModelRM; % if RM to be computed, provide r0 lattice without errors
end

traj=0; % initial correction uses model lattice
% if needed will switch to traj=1 during the loop.
if printouttext
    if traj==1
        disp('correcting trajectory');
    elseif traj==0
        disp('correcting trajectory with model rm')
    end
end

% load  response matrices or compute them
if correctflags(1) % dpp correction
    rmsel=[4 5 6];
else
    rmsel=[4 5];
end

if iscell(ModelRM)
    if printouttext
        disp('computing model trajectory RM')
    end
    
    ModelRM=getresponsematrices(...
        r0,...          %1 AT lattice
        indBPM,...      %2 bpm indexes in at lattice
        indHCor,...     %3 h cor indexes
        indVCor,...     %4 v cor indexes
        [],...     %5 skew cor indexes
        [],...     %6 quad cor indexes
        [],...
        inCOD,...       %7 initial coordinates
        rmsel...        %8 specifiy rm to be computed
        );
    
    if ~correctflags(1) % dpp correction
        ModelRM.TrajHDPP{1}=[];
        ModelRM.TrajVDPP{3}=[];
    end
    
    if printouttext
        disp('computed model trajectory RM')
    end
    
end



RMH=ModelRM.TrajHCor{1};
RMV=ModelRM.TrajVCor{3};
RMHd=ModelRM.TrajHDPP;
RMVd=ModelRM.TrajVDPP;


rclosed=[];

nbpmuprev=0;
nbpmu=0; % bpm used in the correction

lim1=lim(1);

countiter=0;
countstack=0;% chack that correction did not get stack

countatbpm=zeros(size(indBPM));
fraceig=2;
Nforcor=maxcornum;


% loop until a closed orbit is not found
while isempty(rclosed)
    if printouttext
        disp('Search closed orbit');
    end
    if countstack>6 % 6 iterations are limit to continue.
        error(['correction of open trjectory looping at ' num2str(nbpmu) ' bpms'])
    end
    
    % check for closed orbit,
    % if found, end loop,
    % if not found but all bpm, close with last 2 H/V steerers
    if (nbpmu==length(indBPM))
        if printouttext
            disp(['Trajectory closure with last 2 correctors']);
        end
        
        ropen=MatchLast2CorForFirstBPM(ropen, inCOD,...
                                       indBPM, indHCor, indVCor);
        
        COD=findorbit6Err(ropen,indBPM,inCOD);
        
        if ~isempty(find(isnan(COD),1))
            if printouttext
                disp('no closed orbit, but all bpm see signal. Closing with last 2 correctors.'); end;
            
            rclosed=MatchLast2CorForFirstBPM(ropen, inCOD,...
                                             indBPM, indHCor, indVCor);
        end
        
        if isempty(find(isnan(COD),1))
            if printouttext
                disp('Found closed orbit'); end;
            rclosed=ropen; % end loop!
        end
        
    else % less then All BPM see signal ( 1 turn not completed)
        
        countiter=countiter+1;
        
        % get best input coordinates (septum tuning)
        %disp('Scan input coordinates')
        %[inCOD,nelpas]=ScanPosAngle(ropen,r0,indBPM,inCOD,51,3e-3,'scanInCOD');
        %disp(['Best Input Coordinates: [' num2str(inCOD*1e3,'%.2d, ') '] mm, ' num2str(nelpas) ' elem'])
        
        % get current orbit
        [t]=findtrajectory6Err(ropen,indBPM,inCOD);
        
        ox=t(1,:)-reforbit(1,:);
        oy=t(3,:)-reforbit(2,:);
        
        % bpm at wich there is a signal (not NaN and less than 10mm)
        usebpm=~isnan(ox) & ~isnan(oy) & ox<lim(1) & ox>-lim(1) & oy<lim(1) & oy>-lim(1);
        
        %  % stop at first bpm giving NAN
        %  % usebpm=1:(find(isnan(ox) || isnan(oy),1,'first')-1);
        
        % restrict to smaller initial region
        usebpm=indBPM(usebpm)<=max(indBPM(usebpm))*1;
        
        if countiter>1
            nbpmuprev=nbpmu;% store previous value
        end
        
        nbpmu=length(find(usebpm));
        
        
        if nbpmuprev==nbpmu % no more bpm then previous correction
            if countstack>1
                lim1=lim1+0.5e-3; % accept less bpm
                if printouttext
                    disp(['no improvement from last step. increasing lim to: ' num2str(lim1)]); end;
                usebpm=~isnan(ox) & ~isnan(oy) & ox<lim1 & ox>-lim1 & oy<lim1 & oy>-lim1;
                
                if lim1>lim+2e-3
                    if printouttext
                        warning('Could not find closed orbit. Trajectory above 6mm. Aborting trajectory correction'); end;
                    warning('Closed orbit NOT FOUND');
                    rclosed=ropen;
                    return
                end
                
                if lim1>lim+1e-3 && numel(find(usebpm))<320
                    if printouttext
                        disp('update optimal injection point'); end;
                    [inCOD]=Scan2x2DinCOD(ropen,inCOD,101,lim1,[]);
                    
                end
                
                
            end
            countstack=countstack+1;
            
            
            % actions if number of corrector is stuck for >3 iterations
            if countstack==3
                traj=0;
                fraceig=5;% 1/5 eigenvectors
            end
            if countstack==4
                traj=0;
                fraceig=10; % 1/10 eigenvectors
            end
            if countstack==5
                traj=1;% recompute trajectory rm
                if printouttext
                    disp('setting correction to measured trajectory response matrix.');end;
                fraceig=7; % 1/7 eigenvectors
            end
            if countstack==6
                traj=1;% recompute trajectory rm
                if printouttext
                    disp('setting correction to measured trajectory response matrix.'); end;
                fraceig=10;% 1/10 eigenvectors
            end
            
        else
            countstack=0;
            
        end
        
        if countstack==10
            error('Could not find closed orbit. Aborting trajectory correction');
        end
        
        nbpmu=length(find(usebpm));
        if nbpmu==0
            error('NO BEAM AT ALL!')
        end
        
        % all available correctors before last bpm reading signal
        usecorH=indHCor<max(indBPM(usebpm));
        usecorV=indVCor<max(indBPM(usebpm));
        
        nhcu=length(find(usecorH));
        nvcu=length(find(usecorV));
        
        countatbpm(nbpmu)=countatbpm(nbpmu)+1;
        
        if nhcu==0 || nvcu==0
            error('NO BEAM AT ALL!')
        end
        
        if printouttext
            disp(['Trajectory correction:'...
                ' nbpms= ' num2str(nbpmu)...
                ' ncor: ' num2str([nhcu nvcu],'%d, ')]);
        end
        
        % limit to last correctors and bpms
        
        lastbpm=find(usebpm==1,1,'last');
        if lastbpm>Nforcor
            usebpm(1:lastbpm-Nforcor)=0;
        end
        
        lastcor=find(usecorH==1,1,'last');
        if lastcor>Nforcor
            usecorH(1:lastcor-Nforcor)=0;
        end
        
        lastcor=find(usecorV==1,1,'last');
        if lastcor>Nforcor
            usecorV(1:lastcor-Nforcor)=0;
        end
        
        if printouttext
            disp('computing ORM for available trajectory');
        end
        
        bpm=indBPM(usebpm);
        
        X=ox(usebpm)';
        Y=oy(usebpm)';
        
        HK=indHCor(usecorH);
        VK=indVCor(usecorV);
        
        neigSteerer=floor(length(HK)/fraceig);
        fracapply=1;
        
        % ----- CORRECTION ----- %
        
        if printouttext
            disp('H PLANE');
        end;
        corh0=getcellstruct(ropen,'PolynomB',HK,1,1);
        
        if traj==1
            
            ModelRM=getresponsematrices(...
                ropen,...          %1 AT lattice
                bpm,...      %2 bpm indexes in at lattice
                HK,...     %3 h cor indexes
                [],...     %4 v cor indexes
                [],...     %5 skew cor indexes
                [],...     %6 quad cor indexes
                [],...
                inCOD,...       %7 initial coordinates
                rmsel...        %8 specifiy rm to be computed
                );
            if ~correctflags(1) % dpp correction
                ModelRM.TrajHDPP{1}=[];
            end
            RespH=ModelRM.TrajHCor{1};
            RespHd=ModelRM.TrajHDPP;
        elseif traj==0
            RespH=RMH(usebpm,usecorH);
            if correctflags(1) % dpp correction
                RespHd=RMHd(usebpm);
            end
        end
        
        % build RMs
        if      correctflags(1) &&  correctflags(2) % dpp and mean0
            FRH=[ [RespH;ones(1,size(RespH,2))] [RespHd';0] ];
        elseif  correctflags(1) && ~correctflags(2) % dpp no mean 0
            FRH=[ RespH RespHd' ];
        elseif ~correctflags(1) &&  correctflags(2) % no dpp mean0
            FRH=[RespH;ones(1,size(RespH,2))];
        elseif ~correctflags(1) && ~correctflags(2) % no dpp no mean0
            FRH=RespH;
        end
        
        % if mean0 correction add 0 to correction vector
        if correctflags(2)
            dch=qemsvd_mod(FRH,[X;0],neigSteerer)*fracapply;
        else
            dch=qemsvd_mod(FRH,X,neigSteerer)*fracapply;
        end
        
        % if dpp correction separate dpp from correctors
        if correctflags(1)
            %    inCOD(5)=inCOD(5)+dch(end);
            dch=dch(1:end-1);
            deltacor=dch(end);
        end
        
        % limit steerers strengths
        if ~isempty(steererlimit)
            dch(abs(dch)>steererlimit(1))=steererlimit(1);
        end
        
        ropen=setcellstruct(ropen,'PolynomB',HK,corh0-dch,1,1);
        %ropen=atsetfieldvalues(ropen,indrfc,'Frequency',f0-alpha*(deltacor)*f0);%
        %alpha not computed since lattice is open.
            
        if printouttext
            disp('correcting available V trajectory');
            disp('V PLANE');
        end;
        
        corv0=getcellstruct(ropen,'PolynomA',VK,1,1);
        
        if traj==1
            
            ModelRM=getresponsematrices(...
                ropen,...          %1 AT lattice
                bpm,...      %2 bpm indexes in at lattice
                [],...     %3 h cor indexes
                VK,...     %4 v cor indexes
                [],...     %5 skew cor indexes
                [],...     %6 quad cor indexes
                [],...
                inCOD,...       %7 initial coordinates
                rmsel...        %8 specifiy rm to be computed
                );
            if ~correctflags(1) % dpp correction
                ModelRM.TrajVDPP{1}=[];
            end
            RespV=ModelRM.TrajVCor{3};
            RespVd=ModelRM.TrajVDPP;
        elseif traj==0
            RespV=RMV(usebpm,usecorH);
            if correctflags(1) % dpp correction
                RespVd=RMVd(usebpm);
            end
        end
        
        % build RMs
        if  correctflags(2) % no dpp mean0
            FRV=[RespV;ones(1,size(RespV,2))];
        elseif ~correctflags(2) % no dpp no mean0
            FRV=RespV;
        end
        
        % if mean0 correction add 0 to correction vector
        if correctflags(2)
            dcv=qemsvd_mod(FRV,[Y;0],neigSteerer)*fracapply;
        else
            dcv=qemsvd_mod(FRV,Y,neigSteerer)*fracapply;
        end
        
        
        % limit steerers strengths
        if ~isempty(steererlimit)
            dcv(abs(dcv)>steererlimit(2))=steererlimit(2);
        end
        ropen=setcellstruct(ropen,'PolynomA',VK,corv0-dcv,1,1);
   
        
        if printouttext
            disp('correcting available V trajectory'); end;
        
        % get current trajectory
        
        [t]=findtrajectory6Err(ropen,indBPM,inCOD);
        
        oxc=t(1,:);
        oyc=t(3,:);
        
        
        if printouttext
            % display results
            disp(['X: ' num2str(std(ox(usebpm))*1e6,'%3.3f') ' -> '...
                num2str(std(oxc(usebpm))*1e6,'%3.3f') ' um'])
            disp(['Y: ' num2str(std(oy(usebpm))*1e6,'%3.3f') ' -> '...
                num2str(std(oyc(usebpm))*1e6,'%3.3f') ' um'])
        end;
        
    end
    
end

if printouttext
    
    [t]=findtrajectory6Err(ropen,indBPM,inCOD);
    
    ox=t(1,:);
    oy=t(3,:);
    
    [t]=findtrajectory6Err(rclosed,indBPM,inCOD);
    
    oxc=t(1,:);
    oyc=t(3,:);
    
    
    % display results
    disp(['X: ' num2str(std(ox)*1e6,'%3.3f') ' -> '...
        num2str(std(oxc)*1e6,'%3.3f') ' um'])
    disp(['Y: ' num2str(std(oy)*1e6,'%3.3f') ' -> '...
        num2str(std(oyc)*1e6,'%3.3f') ' um'])
end;

% update initial cloed orbit guess
[inCOD]=findorbit6Err(rclosed,1,inCOD);


end