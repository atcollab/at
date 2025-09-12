function deltamax = momentum_aperture_at( ...
                                        THERING, ...
                                        deltalimit, ...
                                        initcoord, ...
                                        delta, ...
                                        precdelta, ...
                                        deltastepsize, ...
                                        splits, ...
                                        split_step_divisor, ...
                                        nturns, ...
                                        varargin ...
                                        )
%momentum_aperture_at recursively offsets the particle energy and checks
%                     for survival over n turns of tracking.
%                     Returns the stable energy boundary.
%
% deltamax ...
%     = momentum_aperture_at( ...
%         THERING,...
%         deltalimit,...       [min max]
%         initcoord,...        [x y] initial coordinate
%         delta,...            current energy offset
%         precdelta,...        previous energy offset
%         deltastepsize,...
%         splits,...           number of times splitting the deltastepsize
%         split_step_divisor,  divides the step size at every split
%         nturns
%         )
%
% ... = momentum_aperture_at(..., 'reforbit',ORBITIN)
%       Use ORBITIN to define the reference. Useful when the closed orbit
%       is not zero.
%
% Adapted routine based on ELEGANT
% 1. Start with delta = 0, i.e., zero momentum offset.
% 2. If the limit has been reached stop, otherwise
%      If the number of step divisions is done, stop. Otherwise ...
%      Track the particle
%      If it survives, increase the energy by one step, and start 2) again.
%      Otherwise, go back one step in energy, and divide the step.
%      Count the number of times the step has been divided.
%      Start 2) with the new step.
%
% Debugging info prints are commented to avoid speed reduction,
%
% The ELEGANT routine:
% https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu53.html
%
% ex: [deltamax]=momentum_aperture_at(THERING,0.1,[10^-6 10^-6],0,0,0.01,3,10,100)
% ex: [deltamin]=momentum_aperture_at(THERING,-0.1,[10^-6 10^-6],0,0,-0.01,3,10,100)

% 2024apr09 oblanco at ALBA : adds reforbit option.

p = inputParser;
addOptional(p,'reforbit',zeros(6,1));
addOptional(p,'verbose',false);
parse(p,varargin{:});
par = p.Results;

verbose = par.verbose;

offset6d = par.reforbit + [initcoord(1); 0; initcoord(2); 0; delta; 0];
if ( delta>=0 && delta<deltalimit) ||  ( delta<=0 && delta>deltalimit )
    if splits>-1
        % track for this delta
        [~, LOSS] =ringpass(THERING,offset6d,nturns);
        if LOSS~=1 % if NOT LOST go to next step
            thedelta = delta+deltastepsize;
            thepreviousdelta = delta;
            thedeltastepsize = deltastepsize;
            thesplits = splits;
        else % if LOST reduce stepsize
             % go back to previous delta center and increase of smaller step
             thedelta = precdelta+deltastepsize/split_step_divisor;
             thepreviousdelta = precdelta;
             thedeltastepsize = deltastepsize/split_step_divisor;
             thesplits = splits - 1;
        end
        if verbose
            if LOSS~=1
                deadalivemsg = 'alive';
            else
                deadalivemsg = 'dead';
            end
            fprintf( ...
                'split%d\tdelta%.5f\tprecdelta%.5f\tdeltastepsize%.5f\t%s\n', ...
                splits, ...
                delta, ...
                precdelta, ...
                deltastepsize, ...
                deadalivemsg ...
                );
        end
        deltamax ...
            = momentum_aperture_at( ...
                THERING, ...
                deltalimit, ... [min max]
                initcoord, ... [x y]
                thedelta, ... % delta center
                thepreviousdelta, ...
                thedeltastepsize, ...
                thesplits, ... % number of splitting
                split_step_divisor, ...
                nturns, ...
                'reforbit',par.reforbit, ...
                'verbose', verbose ...
                );
    else
        % no splitting steps remain
        deltamax=delta-deltastepsize;
    end
else
    % limit reached
    deltamax=delta;
end
return;
