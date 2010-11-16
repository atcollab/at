function [nux,nuz]=nuampls(ring,ampl,xz)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=nuampls(RING,AMPLITUDE)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,3)
%
%	Computes tunes for the specified vertical amplitudes
% e.g.  atnuampl(esrf,0:.0002:0.01)
%
% setting the argument plt to 1 will plot the resulting tunes vs.
% amplitude.

if nargin < 3, xz=1; end
if nargin < 4, plt=0; end
if nargin < 5, pltype='0-'; end
siza=size(ampl);
nampl=prod(siza);
p0=repmat(0.00001*[1;0;1;0;0;0], 1,nampl);
p0(xz,:)=p0(xz,:)+ampl(:)';
p1=ringpass(ring,p0,128);
x1=reshape(p1(1,:)-i*p1(2,:),nampl,128)';
z1=reshape(p1(3,:)-i*p1(4,:),nampl,128)';

nux=reshape(findtune(x1,3),siza);

nuz=reshape(findtune(z1,3),siza);
