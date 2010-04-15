function [nux,nuz]=atnuampl(ring,ampl,xz,plt)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE)
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
siza=size(ampl);
nampl=prod(siza);
p0=repmat(0.00001*[1;0;7;0;0;0], 1,nampl);
p0(xz,:)=p0(xz,:)+ampl(:)';
p1=ringpass(ring,p0,128);
x1=reshape(p1(1,:)-i*p1(2,:),nampl,128)';
z1=reshape(p1(3,:)-i*p1(4,:),nampl,128)';
%nux=reshape(findtune(reshape(p1(1,:),nampl,[])',3),siza);
nux=reshape(findtune(x1,3),siza);
%nux=reshape(findtune(reshape(p1(1,:),nampl,[])'),siza);
%nuz=reshape(findtune(reshape(p1(3,:),nampl,[])',3),siza);
nuz=reshape(findtune(z1,3),siza);
%plot((ampl.*ampl)',[nux-nux(1);nuz-nuz(1)]','o-');
if (plt==1)
    plot((ampl),[nux;nuz],'o-');
    legend('\nu_x','\nu_z');
    grid on
end