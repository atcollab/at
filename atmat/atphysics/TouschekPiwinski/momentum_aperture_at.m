function [deltamax]=momentum_aperture_at(THERING,deltalimit,initcoord,delta,precdelta,deltastepsize,splits,split_step_divisor,nturns)
% function [deltamin, deltamax...
%     ]=momentum_aperture_at(THERING,...
%     deltalimit,... [min max]
%     initcoord,... [x y] initial coordinate
%     delta,...
%     precdelta,...
%     deltastepsize,...
%     splits,... % number of splitting
%     split_step_divisor)  % divide the step size at every split
%
% following the ELEGANT routine:
% Start with ? = 0, i.e., zero momentum offset.
% 2. Track a particle to see if it gets lost. If so, proceed to step 4.
% 3. Increase ? by step size ?? and return to step 2.
% 4. If no splitting steps remain, proceed to the next step. Otherwise:
% (a) Change ? to deltas ? sb??., where ?s is the largest ? for which the particle survived, and sb is the steps_back parameter.
% (b) Divide the step size by split_step_divisor to get a new step size ??.
% (c) Set?=?+??.
% (d) Decrement the ?splits remaining? counter by 1.
% (e) Continue from step 2.
% 5. Stop. The momentum aperture is ?s
%
% ex: [deltamax]=momentum_aperture_at(THERING,0.1,[10^-6 10^-6],0,0,0.01,3,10,100)
% ex: [deltamin]=momentum_aperture_at(THERING,-0.1,[10^-6 10^-6],0,0,-0.01,3,10,100)

%disp([delta splits])

if ( delta>=0 && delta<deltalimit) ||  ( delta<=0 && delta>deltalimit)
        
    if splits>-1
                
        % track for this delta
        
        [~, LOSS] =ringpass(THERING,[initcoord(1) 0 initcoord(2) 0 delta 0]',nturns);
        
        if LOSS~=1 % if NOT LOST go to next step
            
            [deltamax...
                ]=momentum_aperture_at(THERING,...
                deltalimit,... [min max]
                initcoord,... [x y]
                delta+deltastepsize,... % delta center
                delta,...
                deltastepsize,...
                splits,... % number of splitting
                split_step_divisor,...
                nturns);
            
        else % if LOST reduce stepsize
            
            [deltamax...
                ]=momentum_aperture_at(THERING,...
                deltalimit,... [min max]
                initcoord,... [x y]
                precdelta+deltastepsize/split_step_divisor,... % go back to previous delta center and increase of smaller step
                precdelta,...
                deltastepsize/split_step_divisor,...
                splits-1,... % number of splitting
                split_step_divisor,...
                nturns);
            
        end
    else
        
        % no splitting steps remain
        deltamax=delta-deltastepsize;
        
    end
    
else
    % limit reached
    deltamax=delta;
end


return;



