function [s,tune] = beam22(t)
%BEAM22 computes the beam matrix from the 1-turn transfer matrix
%%
%BEAM=BEAM22(T)
% T:    1-turn transfer matrix
% BEAM: envelope matrix
%
%[BEAM,TUNE]=BEAM22(T)
% also returns the tune

cosmu=0.5*trace(t);
sinmu=sign(t(1,2))*sqrt(1-cosmu*cosmu);
alpha=0.5*(t(1,1)-t(2,2))/sinmu;
beta=t(1,2)/sinmu;
s=[beta -alpha; -alpha (alpha*alpha+1)/beta];
tune=atan2(sinmu,cosmu)/2/pi;
end

