function [x,y] =atmakeXYProjectionEllipse( sigxx,sigyy,sigxy )
%ATMAKEXYPROJECTIONELLIPSE Gives points to plot the contour ellipses
%given the second spatial moments <x^2>, <y^2>, and <xy> 

%we solve for effective Twiss parameters and emittance.
%Note that since this is x-y space, this is just a computational tool
%This emittance is not constant around ring.

%first compute epsilon
eps=sqrt(sigxx.*sigyy-sigxy.^2);
%now compute beta,gamma,alpha
beta=sigxx./eps;
gamma=sigyy./eps;
alpha=-sigxy./eps;

t=linspace(0,2*pi);
xh=sqrt(eps)*cos(t);
yh=sqrt(eps)*sin(t);

x=sqrt(beta)*xh;
y=(-alpha/sqrt(beta))*xh+(1/sqrt(beta))*yh;

