function [map_l,map_h]=MomAperture_allRing(THERING,points,varargin)
%MomAperture_allRing returns positive and negative momentum aperture
%                    boundaries where the particle is still alive.
%
%                    The boundary width (i.e. the uncertainty) is equal to
%                    energystep / (splitdivisor^ntimessplit),
%                    meaning that one more step of this size makes the
%                    particle unstable.
%
% [map_l,map_h] ...
%    = MomAperture_allRing(
%                          THERING, ...
%                          POINTS ...
%                         )
% [...] = MomAperture_allRing(..., NTURNS)
%         Tracks over NTURNS to get the momentum aperture. Default 100
%         e.g [dppM,dppP]=MomAperture_allRing(THERING,positions,nturns)
%
% [...] = MomAperture_allRing(..., 'reforbit',ORBITIN)
%         The initial particle coordinates are taken from ORBITIN.
%         Default zeros(6,length(POINTS))
%
% [...] = MomAperture_allRing(..., 'xyinitoffsets',[x y])
%         The transverse offsets to add to the reference orbit.
%         Default 1e-5*ones(length(POINTS),2)
%
% [...] = MomAperture_allRing(..., 'deltalimits',[deltapos deltaneg])
%         The energy offset limits. Default [0.1 -0.1]
%
% [...] = MomAperture_allRing(..., 'initialguess',[posguess negguess])
%         The starting point of the recursive energy offsets.
%         Default [0.0 0.0]
%
% [...] = MomAperture_allRing(..., 'energysteps',[posstep negstep])
%         The positive and negative initial energy steps.
%         Default [0.01 -0.01]
%
% [...] = MomAperture_allRing(..., 'ntimessplit',NSPLIT)
%         The number of recursive calls reducing the step size. Default 2
%
% [...] = MomAperture_allRing(..., 'splitdivisor',SPLITTDIVISOR)
%         The step divisor every time we split the step. Default 10.
%
% [...] = MomAperture_allRing(..., 'verbose',VERBOSE)
%         Print info the current position. Default 0.
%         If set to 1 it will print info at every reference point.
%         If set to 2 it will print info at each energy step.
%
% ex: [map_l,map_h] = MomAperture_allRing(THERING,points,nturns);
% ex: [map_l,map_h] = ...
%   MomAperture_allRing(THERING,points,nturns,'reforbit',findorbit6(THERING,points));

% 2024apr09 oblancog at ALBA : adds parameters, and doc info
map_h=zeros(length(points),1);
map_l=zeros(length(points),1);

% nturns is historically a positional argument, we work around to find it
if mod(numel(varargin),2) == 1
    nturns=varargin{1};
    varargin(1) = [];
else
    nturns=100;
end

p = inputParser;
addOptional(p,'reforbit',zeros(6,length(points)));
addOptional(p,'xyinitoffsets',10e-6*ones(length(points),2));
addOptional(p,'deltalimits', [0.1 -0.1]);
addOptional(p,'initialguess', [0.0 0.0]);
addOptional(p,'energysteps', [0.01 -0.01]);
addOptional(p,'ntimessplit',2);
addOptional(p,'splitdivisor',10);
addOptional(p,'verbose',0);
parse(p,varargin{:});
par = p.Results;

verbose = par.verbose;
deepverbose = false;
if verbose == 2
    deepverbose = true;
end

if ~isequal(size(par.reforbit),[6,length(points)])
    fprintf('ERROR: The size of the reference orbit needs to be (6,number of points)\n');
    return;
end

lengthpoints = length(points);
maps = zeros(lengthpoints,2);
deltalimits = par.deltalimits;
esteps = par.energysteps;
for i=1:length(points)
    %cycle ring
    THERING_cycl=[THERING(points(i):end); THERING(1:points(i)-1)]';

    for j = 1:length(deltalimits)
        % try catch here is used to return a valid value even if the momap
        % calculation has problems. At the end it is very expensive to
        % throw an error and exit, better return zero if everything fails
        try
            maps(i,j)=momentum_aperture_at(...
                        THERING_cycl, ...
                        deltalimits(j), ...
                        par.xyinitoffsets(i,:), ...
                        par.initialguess(j), ...
                        par.initialguess(j), ...
                        esteps(j), ...
                        par.ntimessplit, ...
                        par.splitdivisor, ...
                        nturns, ...
                        'reforbit',par.reforbit(:,i), ...
                        'verbose', deepverbose ...
                        );
        catch
            maps(i,j)=0;
        end
    end
    if verbose
      fprintf('%d turns, point %d of %d, %.0f%%\tdppP:%6.3f\tdppN:%6.3f\n', ...
                  nturns, ...
                  i, ...
                  lengthpoints, ...
                  i/lengthpoints*100, ...
                  maps(i,1), ...
                  maps(i,2) ...
                  );
    end
end
map_l = maps(:,2);
map_h = maps(:,1);
return
