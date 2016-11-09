function re=setTiltAbout(r,pos,coord,PSI)
%
%
% update
%
%see also: SETSHIFT_THERING settilt_THERING_Dipole

numelems = length(pos);

if (numelems ~= length(PSI))
    error('ELEMINDEX, D must have the same number of elements');
end

b_s=1; % speed of electrons!

additive=0;

switch coord
    case 'y_MAD8'
        for i = 1:length(pos)
            
            C = cos(PSI);
            S = sin(PSI);
            T = tan(PSI);
            
            RM = diag([1 1  C(i) 1/C(i) 1  1 ]);
            RM(4,5) = -1/b_s*T(i);
            RM(1,6) = 1/b_s*S(i);
            
            r{pos(i)}.R1 = RM;
            r{pos(i)}.R2 = RM';
            r{pos(i)}.RotAboutY =   PSI(i);
         %   r{pos(i)}.Pitch =   PSI(i);
            
        end
    case 'x_MAD8'
        for i = 1:length(pos)% angle is changing!
            C = cos(PSI);
            S = sin(PSI);
            T = tan(PSI);
            
            RM = diag([ C(i) 1/C(i) 1 1 1  1 ]);
            RM(2,5) = -1/b_s*T(i);
            RM(1,6) = 1/b_s*S(i);
            
            r{pos(i)}.R1 = RM;
            r{pos(i)}.R2 = RM';
         
             r{pos(i)}.RotAboutX =   PSI(i);
       %     r{pos(i)}.Yaw =   PSI(i);
            
        end
        
    case 'y'
        for i = 1:length(pos)
            if additive
                r{pos(i)}.T1(1,2) =  r{pos(i)}.T1(2)  -PSI(i);
                r{pos(i)}.T2(1,2) = r{pos(i)}.T2(2)   +PSI(i);
                r{pos(i)}.T1(1,1) =  r{pos(i)}.T1(1)  +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.T2(1,1) =   r{pos(i)}.T2(1) +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.RotAboutY = r{pos(i)}.RotAboutY +  PSI(i);
           %     r{pos(i)}.Pitch = r{pos(i)}.Pitch +  PSI(i);
                %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
            else
                r{pos(i)}.T1(1,2) =   -PSI(i);
                r{pos(i)}.T2(1,2) =   +PSI(i);
                r{pos(i)}.T1(1,1) =   +r{pos(i)}.Length/2*tan(PSI(i));
                r{pos(i)}.T2(1,1) =   +r{pos(i)}.Length/2*tan(PSI(i));
                r{pos(i)}.RotAboutY =   PSI(i);
             %   r{pos(i)}.Pitch =   PSI(i);
                %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
            end
            
            % % reset error
% r=atsetfieldvalues(r,ind,'T1',[0 0 0 0 0 0]);
% r=atsetfieldvalues(r,ind,'T2',[0 0 0 0 0 0]);
% 
% % rotate
% r=setTiltAbout(r,ind,'x',rotval*ones(size(ind)));
% 
% % displace
% t1=atgetfieldvalues(r,ind,'T1',{1,1});
% t2=atgetfieldvalues(r,ind,'T2',{1,1});
% 
% r=atsetfieldvalues(r,ind,'T1',{1,1},t1-dispval/cos(rotval));
% r=atsetfieldvalues(r,ind,'T2',{1,1},t2+dispval/cos(rotval));

            
        end
    case 'x'
        for i = 1:length(pos)
            if additive
                r{pos(i)}.T1(1,4) = r{pos(i)}.T1(4)  -PSI(i);
                r{pos(i)}.T2(1,4) = r{pos(i)}.T2(4)  +PSI(i);
                r{pos(i)}.T1(1,3) = r{pos(i)}.T1(3)  +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.T2(1,3) = r{pos(i)}.T2(3)  +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.RotAboutX = r{pos(i)}.RotAboutX+  PSI(i);
          %      r{pos(i)}.Yaw = r{pos(i)}.Yaw+  PSI(i);
                %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
            else
                r{pos(i)}.T1(1,4) =  -PSI(i);
                r{pos(i)}.T2(1,4) =  +PSI(i);
                r{pos(i)}.T1(1,3) =  +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.T2(1,3) =  +r{pos(i)}.Length/2*sin(PSI(i));
                r{pos(i)}.RotAboutX =   PSI(i);
         %       r{pos(i)}.Yaw =   PSI(i);
                %DS =+L/2(1-cos(psi) both in and out, ignored. (if including, do not change sign!)
            end
        end
    case 's'
         
        % assign tilt field for easy recovery of the set error values
        %r(pos)=cellfun(@(el,rot)setfield(el,'Tilt',rot),r(pos),num2cell(PSI),'un',0);
        
        r=atsettiltdipole(r,pos,PSI);
    otherwise
        disp('tilt about x, y or s');
end

re=r;
