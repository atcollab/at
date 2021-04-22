function [elemdata0,ringdata,elemdata] = atlinopt6(ring, varargin)
%ATLINOPT6 Performs linear analysis of the lattice
%
% LinData = ATLINOPT6(RING,REFPTS) is a MATLAB structure array with fields
%
%   R1          - 6x6 lattice functions for x motion
%   R2          - 6x6 lattice functions for z motion
%   R3          - 6x6 lattice functions for s motion
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
[a0,vps]=amat(m66);
tunes=mod(angle(vps)/2/pi,1);
damping_rates=-log(abs(vps));

[r0,mu,r1,r2,r3]=r_analysis(a0, ms);

[alpha0,beta0,disp0]=output(r0{:});
[alpha,beta,dispersion]=cellfun(@output,r1,r2,r3,'UniformOutput',false);

elemdata0=struct(...
    'R1',r0{1},'R2',r0{2},'R3',r0{3},...
    'alpha', alpha0,'beta', beta0,...
    'Dispersion', disp0,...
    'mu', zeros(1,3));

ringdata=struct('tunes',tunes,'damping_rates',damping_rates);

elemdata=struct(...
    'R1',r1,'R2',r2,'R3',r3,...
    'alpha', alpha,'beta', beta,...
    'Dispersion', dispersion,...
    'mu', num2cell(unwrap(cat(1,mu{:})),2));

    function [alpha,beta,dispersion]=output(r1,r2,r3)
        % Extract output parameters from R matrices
        alpha=[r1(2,1) r2(4,3)];
        beta=[r1(1,1) r2(3,3)];
        dispersion=r3(1:4,5)/r3(5,5);
    end

    function up = unwrap(p)
        % Unwrap negative jumps in betatron
        jumps = diff([zeros(1,size(p,2));p],1,1) < -1.e-3;
        up = p+cumsum(jumps)*2*pi;
    end

end

