function newring = ataddmpoleerrors(ring,type,newindex,strength,radius,randflag)
%ataddrandmpole adds a random multipole component to all elements of type
%'type' where type can be 'dipole', 'quadrupole', or 'sextupole'
%
%[newring] = ATRANDMPOLE(ring,type,newindex,strength,radius)
%
%ring = input ring
%type = 'dipole', 'quadrupole' or 'sextupole'
%newindex: index of Multipole to add
%strength: strength of the multipole component at the given radius
%radius: reference radius for defining the absolute strength
%if randflag is set to 1, then the errors will be random, Gaussian
%distributed
% The formula for the added errors is
% B^(N)_(n) = radius^(n-N)*b_n/b_N
% It represents the relative field error to the design value at the ref.
% radius
%For example, to add a random octupole error of 1e-4 at 25 mm, relative to all
%quadrupoles:
% newring =ataddrandmpole(ring,'quadrupole',4,.1e-4,.025,1);
%
%See also: ataddmpolecomppoly attiltelem atshiftelem

% first find all elements of the given type.
% then run through them, and 
% create a new element with new polynom, using ataddmpolecomppoly to make the new
% PolyNom with a random strength scaled by strength and
% atelem to make the element.  Now replace the old element with the new.

if (strcmp(type,'dipole'))
    elemindex0=finddipoles(ring);
    elemindex=find(elemindex0);
    refindex=1;
else
    if (strcmp(type,'quadrupole'))
    elemindex0=findquadrupoles(ring);
    elemindex=find(elemindex0);
    refindex=2;
    else
        elemindex=[];
    end
end

newring=ring;
for j=1:length(elemindex)
elmnt=ring{elemindex(j)};
polyB = elmnt.PolynomB;

if(randflag)
    strength = strength*randn;
end
polyB2 = ataddmpolecomppoly(polyB,refindex,newindex,strength,radius);
elmnt.PolynomB=polyB2;
elmnt.MaxOrder=length(polyB2);
newring{elemindex(j)}=elmnt;
end

function quads = findquadrupoles(ring)
dipoles = finddipoles(ring);
isquadrupole=@(elem,polyb) length(polyb) >= 2 && polyb(2)~=0;
quads=atgetcells(ring,'PolynomB',isquadrupole) & ~dipoles;

function dipoles = finddipoles(ring)
isdipole=@(elem,bangle) bangle~=0;
dipoles=atgetcells(ring,'BendingAngle',isdipole);
