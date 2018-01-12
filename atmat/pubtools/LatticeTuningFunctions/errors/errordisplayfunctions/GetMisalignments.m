function [X,Y,T]=GetMisalignments(THERING,varargin)
% this function retrives 3 vectors, for x and y misalignments and tilts
% the vectors are length of THERING

if numel(varargin)==1
    indxerrors=varargin{1};
else
    indxerrors=1:length(THERING);
end

X=zeros(size(indxerrors));
Y=zeros(size(indxerrors));
t1=findcells(THERING(indxerrors),'T1');
if ~isempty(t1)
    X(t1)=-atgetfieldvalues(THERING,indxerrors(t1),'T1',{1,1});
    Y(t1)=-atgetfieldvalues(THERING,indxerrors(t1),'T1',{3,1});
end

T=zeros(size(indxerrors));
%r1=findcells(THERING(indxerrors),'R1');
tiltedelem=find(atgetcells(THERING,indxerrors,'Tilt'))';
rotelem=find(atgetcells(THERING,indxerrors,'RotAboutS'))';

if ~isempty(tiltedelem) %|| ~isempty(rotelem)
  
    T(tiltedelem)=atgetfieldvalues(THERING,indxerrors(tiltedelem),'Tilt');
    T(rotelem)=atgetfieldvalues(THERING,indxerrors(rotelem),'RotAboutS');
    T(T<1e-7 & T>-1e-7)=0;
end
