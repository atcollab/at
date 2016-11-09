function ringerr=SetESRFAlgeAlignmentError(ring,algeexcelfile,filenameout,nseeds,fracerr,plotdataALGE)
%function SetESRFAlgeAlignmentError(...
%     ring,...          lattice AT
%     algeexcelfile,... excel file with sheets: Coords, dR, dL, dZ
%     filenameout,...   output filename = ['ALGEErrorsLattice_' num2str(ierr) '_' filenameout]
%     nseeds,...        vector of seeds to use
%     plotdataALGE)     % flag to plot 3D alignment data
%
% sets alignment errors defined in excel file algeexcelfile. BPM move with
% the lattice
%
%see also: setDSerr setshift_THERING
if nargin<5
fracerr=1;
end

if nargin<6
    plotdataALGE=0;
end

globalpos = xlsread(algeexcelfile,'Coords');
x=globalpos(:,1);
y=globalpos(:,2);
z=globalpos(:,3);

ringlength=findspos(ring,length(ring)+1);

% get s along lattice for later interpolation and plot
s=[cumsum(sqrt(diff(x).^2+diff(y).^2+diff(z).^2)); ringlength];

dR = xlsread(algeexcelfile,'dR')*1e-3; % [m]
dL = xlsread(algeexcelfile,'dL')*1e-3;
dZ = xlsread(algeexcelfile,'dZ')*1e-3;

% close loop
s=[0;s];
dR=[dR(end,:);dR];
dL=[dL(end,:);dL];
dZ=[dZ(end,:);dZ];

% lattice positions
slat=findspos(ring,1:length(ring));
% shift by 4 cells to go to injection cell

celllength=ringlength/32;

slatshift=slat-celllength*4;

[~,indslatsort]=sort(slatshift);

dxlat=interp1(s,dR,slat(indslatsort));
dzlat=interp1(s,dZ,slat(indslatsort));
dslat=interp1(s,dL,slat(indslatsort));

% remove non significant values
dslat(dslat<0)=0;
dslat(isnan(dslat))=0;
dslat(isinf(dslat))=0;
dxlat(isnan(dxlat))=0;
dxlat(isinf(dxlat))=0;
dzlat(isnan(dzlat))=0;
dzlat(isinf(dzlat))=0;

dslat(end-2:end,1)=dslat(2);
dzlat(end-2:end,1)=dzlat(2);
dxlat(end-2:end,1)=dxlat(2);

if nargin<4
    nseeds=1;%:size(dR,2);
end

if find(nseeds>size(dR,2))
    nseeds=size(dR,2);
end

indbpm=unique([find(atgetcells(ring,'Class','Monitor'))' find(atgetcells(ring,'BetaCode','PU'))']);

%for 
ierr=nseeds;
    
    ringerr=ring;
    
%     % bpm offsets  %%% removed this part and modified ApplyBPMErr.m to
%     consider also misalignment errors T1, R1, as offset and rotation.
%     This implements somehow the concept of positioning error of the BPM
%     as different form the offset error. 
%
%     ofx0=atgetfieldvalues(ringerr,indbpm,'Offset',{1,1});
%     if isnan(ofx0), ofx0=zeros(size(indbpm)); end;
%     
%     ofy0=atgetfieldvalues(ringerr,indbpm,'Offset',{1,2});
%     if isnan(ofy0), ofy0=zeros(size(indbpm)); end;
%     
%     ringerr=setcellstruct(ringerr,'Offset',indbpm,ofx0(:)-dxlat(indbpm,ierr),1,1);
%     ringerr=setcellstruct(ringerr,'Offset',indbpm,ofy0(:)-dzlat(indbpm,ierr),1,2);
%     
    % % longitudinal erros
    % ringerr=setDSerr(ringerr,1:length(ring),dslat(:,ierr)'); % since after this comand there is an atsetshift, the T2 will be ignored!
    
    % transverse errors
    ringerr=atsetshift(ringerr,1:length(ring),dxlat(:,ierr)*fracerr,dzlat(:,ierr)*fracerr); %#ok<NASGU>
    
    if ~isempty(filenameout)
        save(['ALGEErrorsLattice_' num2str(ierr) '_' filenameout],'ringerr');
    end
%end




if plotdataALGE
    
    figure;plot3(x,y,z);
    zlabel('z [m]');xlabel('x [m]');ylabel('y [m]');
    saveas(gca,['globalpos3d' filenameout '.fig']);
    saveas(gca,['globalpos3d' filenameout '.jpg']);
    
    figure;plot(s,dR,'-');xlabel('s [m]');ylabel('dR [m]')
    saveas(gca,['dRdeviation' filenameout '.fig']);
    saveas(gca,['dRdeviation' filenameout '.jpg']);
    figure;plot(s,dZ,'-');xlabel('s [m]');ylabel('dZ [m]')
    saveas(gca,['dZdeviation' filenameout '.fig']);
    saveas(gca,['dZdeviation' filenameout '.jpg']);
    figure;plot(s,dL,'-');xlabel('s [m]');ylabel('dL [m]')
    saveas(gca,['dLdeviation' filenameout '.fig']);
    saveas(gca,['dLdeviation' filenameout '.jpg']);
    
    [~,r,z]=cart2pol(x,y,z);
    dR_=dR(2:end,:);
    dL_=dL(2:end,:);
    dZ_=dZ(2:end,:);
    
    [dX,dY]=pol2cart(dR_,dL_.*repmat(r,1,50));
    
    
    figure;
    plot3(repmat(x,1,50)+dX,repmat(y,1,50)+dY,repmat(0*z,1,50)+dZ_,'-');
    hold on;
    plot3(x,y,z,'r-','LineWidth',3);
    zlabel('z [m]');xlabel('x [m]');ylabel('y [m]');
    saveas(gca,['globalpos3dAndErrors' filenameout '.fig']);
    saveas(gca,['globalpos3dAndErrors' filenameout '.jpg']);
    
    
    figure;
    plot3(dR',dL',dZ','.');
    xlabel('dR');ylabel('dL');zlabel('dZ')
    saveas(gca,['errorclouds3D' filenameout '.fig']);
    saveas(gca,['errorclouds3D' filenameout '.jpg']);
    
end


return