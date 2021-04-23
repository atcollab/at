function [elemdata,ringdata] = atlinopt6(ring, varargin)
%ATLINOPT6 Performs linear analysis of the lattice
%
% [ELEMDATA,RINGDATA] = ATLINOPT6(RING,REFPTS)
%
%For circular machines, ATLINOPT6 analyses
%the 4x4 1-turn transfer matrix if radiation is OFF, or
%the 6x6 1-turn transfer matrix if radiation is ON.
%
%For a transfer line, The "twiss_in" intput must contain either:
% - a field 'R', as provided by ATLINOPT6, or
% - the fields 'beta' and 'alpha', as provided by ATLINOPT and ATLINOPT6
%
%ELEMDATA is a structure array with fields:
%   R           - DxDx(D/2) R-matrices
%   M           - DxD transfer matrix M from the beginning of RING
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   mu          - [mux, muy] 	betatron phase advances
%   beta        - [betax, betay]                 1x2 beta vector
%   alpha       - [alphax, alphay]               1x2 alpha vector
%
%RINGDATA is a structure array with fields:
%   TUNES
%   DAMPING_RATES
%
% [...] = ATLINOPT6(...,'orbit',ORBITIN)
%   Do not search for closed orbit. Instead ORBITIN,a 6x1 vector
%   of initial conditions is used: [x0; px0; y0; py0; DP; 0].
%   The sixth component is ignored.
%   This syntax is useful to specify the entrance orbit if RING is not a
%   ring or to avoid recomputing the closed orbit if is already known.
%
% [...] = ATLINOPT6(...,'twiss_in',TWISSIN)
%   Computes the optics for a transfer line.
%       TWISSIN is a scalar structure with fields:
%           R       4x4x2 R-matrices
%                   or
%           beta    [betax0, betay0] vector
%           alpha	[alphax0, alphay0] vector
%
%  See also atx atlinopt

[twiss_in,varargs]=getoption(varargin,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');
[dp,varargs]=getoption(varargs,'dp',[]);
%[ct,varargs]=getoption(varargs,'ct',[]);
[refpts,~]=getargs(varargs,1);

if isempty(twiss_in)        % Circular machine
    [mxx,ms]=build_1turn_map(ring,refpts,check_radiation(ring));
    nv=size(mxx,1);
    dms=nv/2;
else                        % Transfer line
    if isempty(orbitin), orbitin=zeros(6,1); end
    sigma=build_sigma(twiss_in);
    nv=size(sigma,1);
    dms=nv/2;
    [~,ms]=build_1turn_map(ring,refpts,dms >= 3);
    mxx=sigma*jmat(dms);
end

if dms >= 3
    output=@output6;
else
    output=@output4;
end

[a0,vps]=amat(mxx);
tunes=mod(angle(vps)/2/pi,1);
damping_rates=-log(abs(vps));

ms=squeeze(num2cell(ms,[1 2]));
[~,mu,ri]=r_analysis(a0, ms);
mu=num2cell(unwrap(cat(1,mu{:})),2);

ringdata=struct('tunes',tunes,'damping_rates',damping_rates);

% elemdata0=output(zeros(1,3),r0{:});
elemdata=cellfun(output,ms,mu,ri);

    function elemdata=output6(ms,mu,ri)
        % Extract output parameters from R matrices
        elemdata=struct('R',ri,'M',ms,...
        'alpha', [ri(1,2,1) ri(2,4,2)],'beta', [ri(1,1,1) ri(2,3,2)],...
        'Dispersion',ri(1:4,5,3)/ri(5,5,3),...
        'mu', mu);
    end

    function elemdata=output4(ms,mu,ri)
        % Extract output parameters from R matrices
        elemdata=struct('R',ri,'M',ms,...
        'alpha', [ri(1,2,1) ri(2,4,2)],'beta', [ri(1,1,1) ri(2,3,2)],...
        'mu', mu);
    end

    function up = unwrap(p)
        % Unwrap negative jumps in betatron
        jumps = diff([zeros(1,size(p,2));p],1,1) < -1.e-3;
        up = p+cumsum(jumps)*2*pi;
    end

    function [mxx,ms]=build_1turn_map(ring,refpts,is6d)
        % Build the initial distribution at entrance of the transfer line
        if is6d
            if ~isempty(dp),warning('AT:linopt','In 6D, "dp" and "ct" are ignored'); end
            [mxx,ms]=findm66(ring,refpts,'orbit',orbitin,'XYStep',XYStep,'DPStep',DPStep);
        else
            if isempty(dp), dp=0; end
            [mxx,ms]=findm44(ring,dp,refpts,'orbit',orbitin,'XYStep',XYStep);
        end
    end

    function sigma=build_sigma(twiss_in)
        % build the sigma matrix = R1 + R2
        if isfield(twiss_in,'R')
            sigma=sum(twiss_in.R,3);
        else
            slices=num2cell(reshape(1:4,2,2),1);
            v=num2cell(cat(1,twiss_in.alpha,twiss_in.beta),1);
            sigma=zeros(4,4);
            cellfun(@sigma22,v,slices,'UniformOutput',false);
        end
        
        function sigma22(ab,slc)
            alpha=ab(1);
            beta=ab(2);
            sigma(slc,slc)=[beta -alpha;-alpha (1+alpha*alpha)/beta];
        end
    end

end

