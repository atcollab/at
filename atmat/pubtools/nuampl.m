function [nux,nuz]=nuampl(ring,ampl,xz,pl)
%ATNUAMPL	computes tune shift with amplitude
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE)
%
%	Computes tunes for the specified horizontal amplitudes
%
%[NUX,NUZ]=ATNUAMPL(RING,AMPLITUDE,3)
%
%	Computes tunes for the specified vertical amplitudes
% e.g.  atnuampl(esrf,0:.0002:0.01)

if nargin < 3, xz=1; end
siza=size(ampl);
nampl=prod(siza);
p0=repmat(0.00005*[1;0;1;0;0;0], 1,nampl);
p0(xz,:)=p0(xz,:)+ampl(:)';
p1=ringpass(ring,p0,128);
x1=reshape(p1(1,:),nampl,128)';
nux=reshape(findtune(reshape(p1(1,:),nampl,[])',3),siza);
nuz=reshape(findtune(reshape(p1(3,:),nampl,[])',3),siza);
%plot((ampl.*ampl)',[nux-nux(1);nuz-nuz(1)]','o-');
if(pl)
    plot((ampl.*ampl)',[nux-nux(1);nuz-nuz(1)]','o-');
    legend('\nu_x','\nu_z');
    grid on
end