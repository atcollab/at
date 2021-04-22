function [elemdata,ringdata] = atlinopt6(ring, varargin)
%ATLINOPT6 Performs linear analysis of the lattice
%
% [ELEMDATA,RINGDATA] = ATLINOPT6(RING,REFPTS)
%
%ELEMDATA is a structure aray with fields:
%   R1          - 6x6 lattice functions for x motion
%   R2          - 6x6 lattice functions for z motion
%   R3          - 6x6 lattice functions for s motion
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   mu          - [mux, muy] 	betatron phase advances
%   beta        - [betax, betay]                 1x2 beta vector
%   alpha       - [alphax, alphay]               1x2 alpha vector
%
%RINGDATA is a structure array with fields:
%   TUNES
%   DAMPING_RATES

[twiss_in,varargs]=getoption(varargin,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');
[refpts,~]=getargs(varargs,1);

if isempty(twiss_in)        % Circular machine
else                        % Transfer line
end

radiation=check_radiation(ring);
if radiation
    [mxx,ms]=findm66(ring,refpts,'orbit',orbitin,'XYStep',XYStep,'DPStep',DPStep);
    output=@output6;
else
    [mxx,ms]=findm44(ring,0.0,refpts,'orbit',orbitin,'XYStep',XYStep);
    output=@output4;
end
nv=size(mxx,1);
dms=nv/2;
[a0,vps]=amat(mxx);
tunes=mod(angle(vps)/2/pi,1);
damping_rates=-log(abs(vps));

[r0,mu,ri{1:dms}]=r_analysis(a0, ms);
mu=num2cell(unwrap(cat(1,mu{:})),2);

ringdata=struct('tunes',tunes,'damping_rates',damping_rates);

% elemdata0=output(zeros(1,3),r0{:});
elemdata=cellfun(output,mu,ri{:});

    function elemdata=output6(mu,r1,r2,r3)
        % Extract output parameters from R matrices
        elemdata=struct('R1',r1,'R2',r2,'R3',r3,...
        'alpha', [r1(2,1) r2(4,3)],'beta', [r1(1,1) r2(3,3)],...
        'Dispersion', r3(1:4,5)/r3(5,5),...
        'mu', mu);
    end

    function elemdata=output4(mu,r1,r2)
        % Extract output parameters from R matrices
        elemdata=struct('R1',r1,'R2',r2,...
        'alpha', [r1(2,1) r2(4,3)],'beta', [r1(1,1) r2(3,3)],...
        'mu', mu);
    end

    function up = unwrap(p)
        % Unwrap negative jumps in betatron
        jumps = diff([zeros(1,size(p,2));p],1,1) < -1.e-3;
        up = p+cumsum(jumps)*2*pi;
    end

end

