function re=setTiltGirderAbout(r,mag_gr,coord,PSI)
%
%
% See also: atsettilt atsettiltdipole atsetshift

%mag_gr=getMagGroupsFromGirderIndex(r);

pos=[mag_gr{:}];

numelems = length(pos);

if (numelems ~= length(PSI))
    error('ELEMINDEX, D must have the same number of elements');
end


switch coord
    case 'x'
        for i = 1:length(pos)
            
            r{pos(i)}.T1(2) =  -PSI(i);
            r{pos(i)}.T2(2) =   PSI(i);
            r{pos(i)}.T1(1) =  -r{pos(i)}.Length/2*sin(PSI(i));
            r{pos(i)}.T2(1) =   r{pos(i)}.Length/2*sin(PSI(i));
            r{pos(i)}.RotAboutX =   PSI(i);
            %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
           
        end
        r=ThetaPhiGirder(r,mag_gr);
        
    case 'y'
        for i = 1:length(pos)% angle is changing!
            
            r{pos(i)}.T1(4) =  -PSI(i);
            r{pos(i)}.T2(4) =   PSI(i);
            r{pos(i)}.T1(3) =  -r{pos(i)}.Length/2*sin(PSI(i));
            r{pos(i)}.T2(3) =   r{pos(i)}.Length/2*sin(PSI(i));
            r{pos(i)}.RotAboutY =   PSI(i);
            %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
            
        end
        r=ThetaPhiGirder(r,mag_gr);
        
    case  's'
        r= atsettiltdipole(r,pos,PSI);
        r= UniformGirderErrors(r,mag_gr);
    otherwise
        disp('tilt about x, y or s');
end



re=r;
