function [dx,dy]=getDispersion(THERING,bpmindx)
% function [dx,dy]=getDispersion(THERING,bpmindx)
% 
% determines dispersion taking two orbits at plus and minus DE=0.0001
% 


DE=0.001;

Op=findorbit4(THERING,DE,bpmindx);
Opx=Op(1,:);
Opy=Op(3,:);

Om=findorbit4(THERING,-DE,bpmindx);
Omx=Om(1,:);
Omy=Om(3,:);

dx=(Opx-Omx)./(2*DE);
dy=(Opy-Omy)./(2*DE);

return