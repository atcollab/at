function [H,Hv]=CurlyH(RING,dp,ind)
%function [H,Hv]=CurlyH(RING,dp,ind)
%
% computes Curly H (dispersion invariant)
% 
% RING :at lattice
% dp  : energy deviation
% ind : reference positions
% 
% output 
% 
% H  : horizontal dispersion invariant 
% Hv : vertical dispersion invariant
%
% 

%TD=twissring(RING,dp,ind);
[TD,~,~]=atlinopt(RING,dp,ind);

 B = cat(1,TD.beta);
 A = cat(1,TD.alpha);
 %[D,Dv] = getDispersion(RING,ind);
 %[Dp,Dpv] = getDispersionAngle(RING,ind);
 DD = cat(2,TD.Dispersion);
 D=DD(1,:);
 Dv=DD(3,:);
 Dp=DD(2,:);
 Dpv=DD(4,:);
 
 G=(1+A(:,1).^2)./B(:,1);
 Gv=(1+A(:,2).^2)./B(:,2);
 
 H=B(:,1).*Dp'.*Dp' + 2*A(:,1).*D'.* Dp' + G.*D'.* D'; 
 Hv=B(:,2).*Dpv'.*Dpv' + 2*A(:,2).*Dv'.* Dpv' + Gv.*Dv'.* Dv';
 