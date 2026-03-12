function [pnotlost, plost] = atfloodfill(ring, varargin)
% [pnotlost, plost] = floodfill(ring)
%
% Finds the 2D acceptance of the ring using Flood Fill.
%
% Flood fill tracks particles from the exterior to the border of the
% acceptance.
% The lost particles are returned in plost.
% The not lost particles are returned in pnotlost.
%
% Parameters:
%   ring:       AT lattice.
% 
% Keyword Arguments:
%   nturns:     Number of turns for the tracking. Default: 1000
%   window:     Min and max coordinate range
%               [Axis1min,Axis1max,Axis2min,Axis2max].
%               Default [-10e-3,10e-3,-10e-3,10e-3].
%               Axis1 and Axis2 are defined by 'axes'.
%   gridsize:   Number of steps per axis. Default [10,10].
%   axes:       Indexes of axes to be scanned. Default is [1,3], i.e. x-y.
%   sixdoffset: Offset to be added. Default zeros(6,1).
%               This is useful to study off-axis acceptance on any plane,
%               or off-momentum acceptance by adding dp to the 5th coord.,
%               to track particles on the closed orbit, or to add
%               a small deviation to the tracked coordinates,
%               e.g. [10e-5 10e-5] in the transverse planes.
%   verbose:    Print extra info. Default 0.
% 
% Returns:
%     pnotlost: array of size (2,n_not_lost) containing the 2D offsets of
%               n_not_lost tracked particles that survived.
%     plost:    array of size (3,n_lost) containing the 2D offsets of
%               n_lost trackedp articles that did not survive; and
%               in the third row the turn on which each particle is lost.
%
% Example:
%     [pnl, pl] = atfloodfill(THERING, nturns=500)
%
%
% Notes:
% This method is recomended for single or low number of CPUs, and,
% it does not scale well for parallel computing.
%
% Based on the article,
%   B. Riemann, M. Aiba, J. Kallestrup, and A. Streun, "Efficient
%   algorithms for dynamic aperture and momentum acceptance
%   calculation in synchrotron light sources", Phys. Rev. Accel.
%   Beams, vol. 27, no. 9, p. 094 002, 2024.
%   doi:10.1103/PhysRevAccelBeams.27.094002

% Author : E. Serra,  UAB and ALBA,  2025 original version in python
%                                    See. IPAC2025, MOPB065.
%                                    APPLICATION OF FAST ALGORITHMS TO
%                                    CALCULATE DYNAMIC AND
%                                    MOMENTUM APERTURE TO THE DESIGN OF
%                                    ALBA II.
%                                    doi: 10.18429/JACoW-IPAC25-MOPB065
% Edited : O. Blanco, ALBA,          2025 matlab version

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse optional arguments
    p = inputParser;
    addOptional(p,'nturns',1000);
    addOptional(p,'window',[-10e-3,10e-3,-10e-3,10e-3]);
    addOptional(p,'gridsize',[10,10]);
    addOptional(p,'axes',[1,3]);
    addOptional(p,'sixdoffset',zeros(6,1));
    addOptional(p,'verbose',0);
    parse(p,varargin{:});
    par = p.Results;

    % Initialize variables
    nturns = par.nturns;
    window = par.window;
    gridsize = par.gridsize;
    sixdoffset = par.sixdoffset;
    verbose = par.verbose;
    axes = par.axes;

    % Initialize output in case we return earlier
    pnotlost = zeros(2,0);
    plost = zeros(3,0);

    if verbose, fprintf('Flood fill starts.\n'); end

    if verbose, fprintf('Max. number of turns: %d\n',nturns); end

    % Create the grid
    nx = gridsize(1);
    ny = gridsize(2);
    npart = nx*ny;
    
    % Check the window dimensions
    if (window(1) == window(2) || window(3) == window(4)) && verbose
        fprintf('Window is quite narrow.\n');
        return;
    end
    % Check the grid
    if any(gridsize < 2)
        fprintf('Horizontal and vertical gridsize is too small.\n');
        return;
    end
    xvals = linspace(min(window(1:2)),max(window(1:2)),nx);
    yvals = linspace(min(window(3:4)),max(window(3:4)),ny);
    points = zeros(2,npart);

    if verbose, fprintf('Number of grid points: %d.\n',npart); end

    % Set the order
    ii = 1;
    for ix = 1:nx
        for iy = 1:ny
            points(:,ii) = [xvals(ix),yvals(iy)]';
            ii = ii + 1;
        end
    end

    % Particle coordinates
    particles = zeros(6,npart);
    particles = particles + sixdoffset;
    particles(axes(1),:) = particles(axes(1),:) + points(1,:);
    particles(axes(2),:) = particles(axes(2),:) + points(2,:);

    % Track the particles for nturns using the flood-fill algorithm.
    % the_queue is replaced by an array with a queue_counter in matlab.
    % The order is fixed to explore the top, left, bottom and right sides.
    the_queue = zeros(npart,1); 
    the_queue(1           :(ny-2)           ) = ((nx-1)*ny+2):(nx*ny-1); %r
    the_queue((ny-1)      :(nx+ny-2)        ) = ((nx-1):-1:0)*ny+1; %b
    the_queue((nx+ny-1)   :(nx+ny+ny-4)     ) = 2:(ny-1); %l
    the_queue((nx+ny+ny-3):(2*nx + 2*(ny-2))) = (nx:-1:1)*ny; %t
    queue_counter = 2*nx + 2*(ny-2);

    % particles that have been tracked
    idxpartdone = false(npart,1);

    % output
    plost = cell(npart,1);
    pnotlost = cell(npart,1);

    if verbose, fprintf('Tracking... '); end
    while queue_counter > 0
        % take one index from the_queue and track that particle
        i = the_queue(queue_counter);
        queue_counter = queue_counter - 1;
        % check if valid and not done
        if (1 <= i && i <= npart) && ~idxpartdone(i)
            idxpartdone(i) = true;
            [~, ~, ~, lossinfo] = ringpass(ring, particles(:,i), nturns);

            if lossinfo.lost
                % if lost, add new points
                newidx = [i+1 i-1 i+ny i-ny];
                maskidx = newidx > 0 & newidx <= npart;
                newidx = newidx(maskidx);
                nextpoints = newidx(~idxpartdone(newidx));
                queue_counter_aux = queue_counter+length(nextpoints);
                the_queue(queue_counter+1:queue_counter_aux) = nextpoints;
                queue_counter = queue_counter_aux;

                % save point
                plost{i} = [ ...
                            particles(axes(1),i); ...
                            particles(axes(2),i); ...
                            lossinfo.turn ...
                           ];
            else
                % if not lost, save point
                pnotlost{i} = [ ...
                                particles(axes(1),i); ...
                                particles(axes(2),i); ...
                              ];
            end
        end
    end

    if verbose, fprintf(' done.\n'); end

    % remove empty cells
    plost = plost(~cellfun('isempty',plost));
    pnotlost = pnotlost(~cellfun('isempty',pnotlost));

    % reshape output
    plost = reshape(cell2mat(plost),3,[]);
    pnotlost = reshape(cell2mat(pnotlost),2,[]);

    if isempty(pnotlost) && verbose
        fprintf('No surviving particles.\n');
    end
end
