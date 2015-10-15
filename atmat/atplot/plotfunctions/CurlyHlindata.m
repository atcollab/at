function [H,Hv]=CurlyHlindata(lindata)
%function [H,Hv]=CurlyHlindata(lindata)
%
% computes Curly H (dispersion invariant)
% 
% lindata is the first ouptut of [lindata,~,~]=atlinopt(...)
% (include at least 2 ouptut arguments for dispersion computation)
% 
% output 
% 
% H  : horizontal dispersion invariant 
% Hv : vertical dispersion invariant
%
% 


 B = cat(1,lindata.beta);
 A = cat(1,lindata.alpha);
 %[D,Dv] = gelindataispersion(RING,ind);
 %[Dp,Dpv] = gelindataispersionAngle(RING,ind);
 DD = cat(2,lindata.Dispersion);
 D=DD(1,:);
 Dv=DD(3,:);
 Dp=DD(2,:);
 Dpv=DD(4,:);
 
 G=(1+A(:,1).^2)./B(:,1);
 Gv=(1+A(:,2).^2)./B(:,2);
 
 H=B(:,1).*Dp'.*Dp' + 2*A(:,1).*D'.* Dp' + G.*D'.* D'; 
 Hv=B(:,2).*Dpv'.*Dpv' + 2*A(:,2).*Dv'.* Dpv' + Gv.*Dv'.* Dv';
