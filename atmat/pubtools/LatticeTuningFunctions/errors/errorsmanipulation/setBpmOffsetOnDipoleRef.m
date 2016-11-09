function rerr=setBpmOffsetOnDipoleRef(rerr)
% % set bpm on curve defined by dipole misalignments
%
%see also: ApplyErrorsRand 

inddip=findcells(rerr,'BendingAngle');
indbpm=findcells(rerr,'Class','Monitor');

[X,Y,~,T,~,~,bpmerr]=GetExistingErrors(rerr,inddip);
    
indbpm=[indbpm(end) indbpm indbpm(1)]; % last and first bpm
    
sposBPM=findspos(rerr,indbpm);
sposDip=findspos(rerr,inddip);

offx=interp1(sposDip,X,sposBPM);offx(isnan(offx))=0;
offy=interp1(sposDip,Y,sposBPM);offy(isnan(offy))=0;
offr=interp1(sposDip,T,sposBPM);offr(isnan(offr))=0;

rerr=setcellstruct(rerr,'Offset',indbpm(2:end-1),bpmerr.offsetx-offx(2:end-1),1,1);
rerr=setcellstruct(rerr,'Offset',indbpm(2:end-1),bpmerr.offsetx-offy(2:end-1),1,2);
rerr=setcellstruct(rerr,'Rotation',indbpm(2:end-1),bpmerr.rotation-offr(2:end-1)*0,1,1);% NO ROTATION!! *0

disp('BPMs on dipole positions curve');

return

