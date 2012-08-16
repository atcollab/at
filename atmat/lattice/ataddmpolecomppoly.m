function [PolynomOut] = ataddmpolecomppoly(Polynom,refindex,newindex,strength,radius,varargin)
%ataddmpolecomppoly adds a multipole component to an existing polynomial,
%scaling it by the reference multipole value and a radius to set the length
%scale
%[PolynomOut] = ATADDMPOLECOMPPOLY(Polynom,index,strength,radius)
%
%Polynom = given field polynomial that we want to add a component to
%refindex = index of magnet: 1 for dipole, 2 for quadrupole
%newindex: index of Multipole to add
%strength: strength of the multipole component at the given radius
%radius: reference radius for defining the absolute strength
% B^(N)_(n) = radius^(n-N)*b_n/b_N
%
%See also: attiltelem, atshiftelem, atmodelem
if (Polynom(refindex)==0)
    if nargin==5
        error('reference polynomial value is 0');
    else
       refvalue=varargin{1};
    end 
  
else
    
   refvalue = Polynom(refindex);
end
    
    
    if(newindex > length(Polynom))
        PolynomOut=[Polynom zeros(1,newindex-length(Polynom))];
    else
        PolynomOut=Polynom;
    end
    
    
    val=power(radius,newindex-refindex)*strength/refvalue;
    PolynomOut(newindex)=PolynomOut(newindex)+val;
end
