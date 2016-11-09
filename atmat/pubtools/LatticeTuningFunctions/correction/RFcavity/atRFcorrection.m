function [...
    rcor,...
    inCOD,...
    fc....
    ]=atRFcorrection(...
    rerr,...
    indBPM,...
    inCOD,...
    neigSteerer,...
    scalefactor,...
    ModelRM,...
    reforbit,...
    printouttext)
% function [...
%    rcor,...           1) corrected lattice
%    inCOD,...          2) initial COD 
%    fc...              3) required frequency change
%    ]=atRFcorrection(...
%     rerr,...          1) AT lattice to correct
%     indBPM,...        2) Nbx1 bpm indexes
%     inCOD,...         3) 6x1 initial COD guess
%     neigSteerer,...   4) 1xNiter eigenvectors for correction dpp (default: [1])
%     scalefactor,...   5) scale factor to correction (default: 1.0)
%     ModelRM,...       6) ModelRM.Orb(H/V)DPP = 6x1 array of orbit
%                          response to dpp
%                          if [] compute RM (default: [])
%     reforbit,...      7) 2xNbpm reference orbit to correct to (default 0*2xNb)
%     printouttext      8) if 1 or true, display rms orbit
%     )
%
% corrects RF frequency and TimeLag to have 0 timeLag at cavities and
% minimize the dispersive contribution to orbit ( final <dpp>_ring =0 )
%
% features impelemented:
% ddp correction
% 6D orbit with BPM errors
% initial coordinate
% correction to reference orbit refx refy
% use atsetfieldvalues, atgetcells
%
%
%see also: qemsvd_mod findorbit6Err getresponsematrices



% response matrix kicks
kval=1e-5;
delta=1e-3;

alpha=mcf(rerr);
indrfc=find(atgetcells(rerr,'Frequency'));

% default arguments
if nargin<8
    printouttext=true;
end

if nargin<2
    if printouttext
        disp('get BPM and Correctors indexes'); end;
    indBPM=finc(atgetcells(rerr,'Class','Monitor'));
end

if nargin<3
    inCOD=[0 0 0 0 0 0]';
end

if nargin<4
    neigSteerer=1;
end


if nargin<5
    if printouttext
        disp(' --- scale set to 1.0'); end;
    scalefactor=1.0;
end

if nargin<6
    if printouttext, disp(' --- computing orbit Response matrix'); end;
    ModelRM=[];
end

if nargin<7
    if printouttext, disp(' --- reference orbit = 0'); end;
    reforbit=zeros(2,length(indBPM));
end

if scalefactor<0 || scalefactor>1
    if printouttext
        disp(' --- scale factor out of range. Set to 1.0'); end;
    scalefactor=1.0;
end

% load or compute response matrix
if isempty(ModelRM)
    % get orbit RM
    if printouttext
        disp('get orbit RM'); end;
    
        ModelRM=getresponsematrices(...
            rerr,...          %1 AT lattice
            indBPM,...      %2 bpm indexes in at lattice
            [],...     %3 h cor indexes
            [],...     %4 v cor indexes
            [],...     %5 skew cor indexes
            [],...     %6 quad cor indexes
            [],...
            inCOD,...       %7 initial coordinates
            [3]...      %8 specifiy rm to be computed
            );
   
end

% kval=ModelRM.kval;
dppH=ModelRM.OrbHDPP;
dppV=ModelRM.OrbVDPP;
% delta=ModelRM.delta;

% get intiial RF
f0=rerr{indrfc(1)}.Frequency;
tlag0=atgetfieldvalues(rerr,indrfc,'TimeLag');
    
% get initial orbit
o=findorbit6Err(rerr,indBPM,inCOD);
ox0=o(1,:);
oy0=o(3,:);

%rerr0=rerr;

% iterate correction
Niter=size(neigSteerer,1);
for iter=1:Niter
    
    if printouttext
        disp(['RF correction iter ' num2str(iter,'%d, ') 'n-eig: ' num2str(neigSteerer(iter,:),'%d, ')]);
    end
    
    % intial RF frequency
    f=rerr{indrfc(1)}.Frequency;
    tlag=atgetfieldvalues(rerr,indrfc,'TimeLag');
        
    % get current orbit
    o=findorbit6Err(rerr,indBPM,inCOD);
    ox=o(1,:);
    
    % subtract reference orbit
    ox=ox-reforbit(1,:);
    
    % build RMs
    
    RMH=dppH';
    
    % compute correction
    dd=qemsvd_mod(RMH,ox',neigSteerer(iter,1));
    
   
    % set RF correction
    rcor=atsetfieldvalues(rerr,indrfc,'Frequency',f-alpha*(-dd)*f);
    
    % reset TimeLag
    orb = findorbit6(rcor,indrfc,inCOD);
    rcor=atsetfieldvalues(rcor,indrfc,'TimeLag',tlag-orb(6,:)');
    
    if printouttext
        disp(['Delta RF : ' num2str(-alpha*(dd)*f) ' Hz']);
        disp(['Delta TimeLag : ' num2str(-orb(6,:)*1e6) ' um']);
    end
    
    % lattice start point for next iteration
    rerr=rcor;
end

% get current orbit
o=findorbit6Err(rcor,indBPM,inCOD);
oxc=o(1,:);

% intial RF frequency
fc=rcor{indrfc(1)}.Frequency;
tlagc=atgetfieldvalues(rcor,indrfc,'TimeLag');
    
inCOD=findorbit6Err(rcor,1,inCOD);

if printouttext
    % display results
    disp(['      before' '    ' '-->' '    ' 'after'])
    disp(['oX: ' num2str(std(ox0-reforbit(1,:))*1e6,'%3.3f') ' -> ' num2str(std(oxc-reforbit(1,:))*1e6,'%3.3f') 'um']);
    disp(['freq: ' num2str(f0,'%3.3f') ' -> ' num2str(fc,'%3.3f') ' Hz ' ' ' num2str(fc-f0,'%3.3f') ' Hz ']);
    disp(['tlag: ' num2str(tlag0','%3.3f, ') ' -> ' num2str(tlagc','%3.3f, ') ' m' ' ' num2str(tlagc'-tlag0','%3.3f') ]);
   disp(['dpp: ' num2str(inCOD(5))])
end
