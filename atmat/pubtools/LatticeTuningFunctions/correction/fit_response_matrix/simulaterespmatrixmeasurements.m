function [rmvec,...
    respvectornorm,...
    respvectorskew...
    ]=simulaterespmatrixmeasurements(...
    r,...
    inCOD,...
    indBPM,...
    indHCor,...
    indVCor,...
    kind,...
    ww,...
    msg)
%[rmvec,
% respvector,
% respvectorskew
% ]=simulaterespmatrixmeasurements(
% r,            1) AT lattice
% inCOD,        2) intial coordinates guess 6x1
% indBPM,       3) bpm indexes
% indHCor,      4) horizontal correctors indexes
% indVCor,      5) vertical correctors indexes     
% kind,         6) rm kind, see below for details
% ww,           7) weigths vector
% msg)          8) message to display after RM computation
%
% prepares response matrix vecotor for error fit
%
% kind may be: Full (ALL RM), Skew (Off Diag), Norm (On Diag). 
%
% features:
% - 6D computation
% - BPM errors considered
% - matrices output in m/rad
% 
%see also:  findrespmat

kval=6.2500e-07;
delta=1e-4;

alpha=mcf(r);
indrfc=find(atgetcells(r,'Frequency'));
f0=r{indrfc(1)}.HarmNumber*PhysConstant.speed_of_light_in_vacuum.value/findspos(r,length(r)+1);            
            
% get tune chrom orbit and dispersion
[l,t,ch]=atlinopt(r,0,indHCor);
bx=arrayfun(@(a)a.beta(1),l);
by=arrayfun(@(a)a.beta(2),l);

o=findorbit6Err(r,indBPM,inCOD);
Ox=o(1,:);
Oy=o(3,:);
% d=finddispersion6Err(r,indBPM,indrfc,alpha,delta,inCOD);
% Dx=d(1,:);
% Dy=d(3,:);
 
Dx=getdisph6D(r,indBPM,indrfc,alpha,delta,inCOD)'; % [m/Hz] *Hz
Dy=getdispv6D(r,indBPM,indrfc,alpha,delta,inCOD)'; % [m/Hz] *Hz

% get bpm resolution value 
bpmresx=atgetfieldvalues(r,indBPM,'Reading',{1,1});
bpmresy=atgetfieldvalues(r,indBPM,'Reading',{1,2});
LH=atgetfieldvalues(r,indHCor,'Length',{1,1});
LV=atgetfieldvalues(r,indVCor,'Length',{1,1});

if nargin<8
    msg='Computed Response matrix vector';
end

if nargin<7
    ww=ones(10,1);
end

if nargin<6
    kind='full';
end

disp(msg);

ormH=findrespmat(r,indBPM,indHCor,kval./sqrt(bx),'PolynomB',1,1,'findorbit6Err',inCOD);
ormV=findrespmat(r,indBPM,indVCor,kval./sqrt(by),'PolynomA',1,1,'findorbit6Err',inCOD);

OH=-ormH{1}./repmat(kval./sqrt(bx).*LH',length(indBPM),1);%./repmat(bpmresx,length(indHCor),1)';
OV=ormV{3}./repmat(kval./sqrt(by).*LV',length(indBPM),1);%./repmat(bpmresy,length(indVCor),1)';
OHV=ormH{3}./repmat(kval./sqrt(bx).*LH',length(indBPM),1);%./repmat(bpmresy,length(indHCor),1)';
OVH=-ormV{1}./repmat(kval./sqrt(by).*LV',length(indBPM),1);%./repmat(bpmresx,length(indVCor),1)';

respvectornorm=[ OH(:) ;... H orm
                 OV(:) ;... V orm
               ]; % response in a column vector

respvectorskew=[ OHV(:) ;... H orm
                 OVH(:) ;... V orm
                ]; % response in a column vector

rm=[OH, OVH;... H orm
    OHV, OV;... V orm
    ];

% select wich rm is the first output.
switch kind
    case {'vdisp','VDISP','VDisp'}
        rmvec=Dy';
    case {'hdisp','HDISP','HDisp'}
        rmvec=Dx';
    case {'full','FULL','Full'}
        rmvec=rm(:);
    case {'fulldisp','FULLDISP','FullDisp'}
        rmvec=[rm(:);ww(1)*Dx';ww(2)*Dy';ww(3)*t'];
    case {'fulldisporb','FULLDISPORB','FullDispOrb'}
        rmvec=[rm(:);ww(1)*Dx';ww(2)*Dy';ww(3)*Ox';ww(4)*Oy';ww(5)*t'];
    case {'fulldisporbchrom','FULLDISPORBCHROM','FullDispOrbChrom'}
        rmvec=[rm(:);ww(1)*Dx';ww(2)*Dy';ww(3)*Ox';ww(4)*Oy';ww(5)*t';ww(6)*ch'];
    case {'skew','SKEW','Skew'}
        rmvec=respvectorskew(:);
    case {'norm','NORM','Norm'}
        rmvec=respvectornorm(:);
    case {'skewdisp','SKEWDISP','SkewDisp'}
        rmvec=[respvectorskew(:);ww(1)*Dy'];
    case {'normdisp','NORMDISP','NormDisp'}
        rmvec=[respvectornorm(:);ww(1)*Dx';ww(2)*t'];
    case {'skewdisporb','SKEWDISPORB','SkewDispOrb'}
        rmvec=[respvectorskew(:);ww(1)*Dy';ww(3)*Ox';ww(3)*Oy'];
    case {'normdisporb','NORMDISPORB','NormDispOrb'}
        rmvec=[respvectornorm(:);ww(1)*Dx';ww(2)*Ox';ww(3)*Oy';ww(4)*t'];
end

%rmvec(isnan(rmvec))=0;
%rmvec(isinf(rmvec))=0;

return

