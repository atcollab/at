function lindata = atlinopt6(ring, varargin)
%ATLINOPT6 Performs linear analysis of the lattice
%
% LinData = ATLINOPT6(RING,REFPTS) is a MATLAB structure array with fields
%
%   B1          - 6x6 lattice functions for x motion
%   B2          - 6x6 lattice functions for z motion
%   B3          - 6x6 lattice functions for s motion
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   mu          - [mux, muy] 	betatron phase advances
%   beta        - [betax, betay]                 1x2 beta vector
%   alpha       - [alphax, alphay]               1x2 alpha vector

check_radiation(ring, true);    % Ensure radiation is ON

[twiss_in,varargs]=getoption(varargin,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');
[refpts,~]=getargs(varargs,1);

if isempty(twiss_in)        % Circular machine
else                        % Transfer line
end

[m66,ms]=findm66(ring,refpts,'orbit',orbitin,'XYStep',XYStep,'DPStep',DPStep);

[mu, b1, b2, b3]=bk_analysis(m66, ms);

[alpha,beta,dispersion,mu]=cellfun(@output,mu,b1,b2,b3,'UniformOutput',false);

lindata=struct(...
    'B1', b1,...
    'B2', b2,...
    'B3', b3,...
    'alpha', alpha,...
    'beta', beta,...
    'Dispersion', dispersion,...
    'mu', num2cell(unwrap(cat(1,mu{:})),2));

    function [alpha,beta,dispersion,mu]=output(mu3,b1,b2,b3)
        % Extract output parameters from Bk matrices
        alpha=[b1(2,1) b2(4,3)];
        beta=[b1(1,1) b2(3,3)];
        dispersion=b3(1:4,5)/b3(5,5);
        mu=mu3(1:3);
    end

    function up = unwrap(p)
        % unwrap negative jumps in betatron
        jumps = diff([zeros(1,size(p,2));p],1,1) < -1.e-3;
        up = p+cumsum(jumps)*2*pi;
    end

end

