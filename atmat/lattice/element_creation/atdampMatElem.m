function elem=atdampMatElem(fname,ring)
%   atdampMatElem creates an element that applies the global damping matrix
%   from the ring
% ring should be without radiation passmethods in the bends


[ringrad,radi,cavi,energy]=atradon(ring);

m66_norad=findm66(ring);
m66_rad=findm66(ringrad);

m66damp=inv(m66_norad)*m66_rad;

elem=atM66(fname,m66damp);