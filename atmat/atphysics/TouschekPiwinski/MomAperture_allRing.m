function [map_l,map_h]=MomAperture_allRing(THERING,points,varargin)
% all Ring momentum aperture
%points=1:10:length(THERING);
map_h=zeros(length(points),1);
map_l=zeros(length(points),1);

if numel(varargin)>0
    nturns=varargin{1};
else
    nturns=100;
end

for i=1:length(points)
   % disp([i length(points) i/length(points)*100])
    %cycle ring
     THERING_cycl=[THERING(points(i):end); THERING(1:points(i)-1)]';
        try
            map_h(i)=momentum_aperture_at(THERING_cycl,+0.1,[10^-6 10^-6],0,0,+0.01,2,10,nturns);
        catch
            map_h(i)=0;
        end
        
        try
            map_l(i)=momentum_aperture_at(THERING_cycl,-0.1,[10^-6 10^-6],0,0,-0.01,2,10,nturns);
        catch
            map_l(i)=0;
        end
       

end


return