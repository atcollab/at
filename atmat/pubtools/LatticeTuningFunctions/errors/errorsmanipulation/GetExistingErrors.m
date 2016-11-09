function [X,Y,S,T,R,P,bpmerrors]=GetExistingErrors(THERING,varargin)
% this function retrives 6 vectors, for x, y,s misalignments,
% Tilts (rot about s), Roll (rot about x) and Pinch (rot about y)
%
% the vectors are length of THERING or of the optional second input
%
% in future also bpm errors and girder errors will be added.
%
%see also: ApplyErrorRand SetExistingErrors

if numel(varargin)==1
    indxerrors=varargin{1};
else
    indxerrors=1:length(THERING);
end

if nargout>=2
    % X and Y misal
    X=zeros(size(indxerrors));
    Y=zeros(size(indxerrors));
    t1=findcells(THERING(indxerrors),'T1');
    if ~isempty(t1)
        X(t1)=-getcellstruct(THERING(indxerrors),'T1',t1,1);
        Y(t1)=-getcellstruct(THERING(indxerrors),'T1',t1,3);
    end
end
if nargout>=3
    % S misal
    S=zeros(size(indxerrors));
    s1=findcells(THERING(indxerrors),'DeltaS');
    if ~isempty(s1)
        S(s1)=-getcellstruct(THERING(indxerrors),'DeltaS',s1,1);
    end
end
if nargout>=4
    % tilt about s-axis
    T=zeros(size(indxerrors));
    tiltedelem=findcells(THERING(indxerrors),'RotAboutS');
    if ~isempty(tiltedelem)
        T(tiltedelem)=getcellstruct(THERING(indxerrors),'RotAboutS',tiltedelem);
        T(T<1e-7 & T>-1e-7)=0;
    end
end
if nargout>=5
    % roll about x-axis
    R=zeros(size(indxerrors));
    tiltedelem=findcells(THERING(indxerrors),'RotAboutX');
    if ~isempty(tiltedelem)
        R(tiltedelem)=getcellstruct(THERING(indxerrors),'RotAboutX',tiltedelem);
        R(R<1e-7 & R>-1e-7)=0;
    end
end
if nargout>=6
    % pinch about y-axis
    P=zeros(size(indxerrors));
    tiltedelem=findcells(THERING(indxerrors),'RotAboutY');
    if ~isempty(tiltedelem)
        P(tiltedelem)=getcellstruct(THERING(indxerrors),'RotAboutY',tiltedelem);
        P(P<1e-7 & P>-1e-7)=0;
    end
end
if nargout>=7
    tiltedelem=findcells(THERING,'Class','Monitor');
    bpmerrors.offsetx=zeros(size(tiltedelem));
    bpmerrors.offsety=zeros(size(tiltedelem));
    bpmerrors.rotation=zeros(size(tiltedelem));
    
    if ~isempty(tiltedelem)
        
try
    ox=cell2mat(cellfun(@(a)a.Offset(1,1),THERING(tiltedelem),'un',0));
        oy=cell2mat(cellfun(@(a)a.Offset(1,2),THERING(tiltedelem),'un',0));
        or=cell2mat(cellfun(@(a)a.Rotation(1,1),THERING(tiltedelem),'un',0));
catch
    ox=zeros(size(tiltedelem));
    oy=zeros(size(tiltedelem));
    or=zeros(size(tiltedelem));
    
end
        ox(isnan(ox))=0;
        oy(isnan(oy))=0;
        or(isnan(or))=0;
        
        bpmerrors.offsetx=ox(:)';
        bpmerrors.offsety=oy(:)';
        bpmerrors.rotation=or(:)';
        
    end
end
return
