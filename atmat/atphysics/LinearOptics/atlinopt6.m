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
%   SPos        - longitudinal position [m]
%   ClosedOrbit - 6x1 closed orbit vector with components
%                 x, px, y, py, dp/d, ct (momentums, NOT angles)
%   R           - DxDx(D/2) R-matrices
%   A           - DxD A-matrix
%   M           - DxD transfer matrix M from the beginning of RING
%   beta        - [betax, betay]                 1x2 beta vector
%   alpha       - [alphax, alphay]               1x2 alpha vector
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   mu          - [mux, muy] 	Betatron phase advances
%   W           - [Wx, Wy]      Chromatic amplitude function [3]
%
%RINGDATA is a structure array with fields:
%   tune            Fractional tunes
%   damping_time    Damping times [s]
%   chromaticity    Chromaticities
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
%  REFERENCES
%   [1] Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
%   [2] Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
%       Published 3 February 2006
%   [3] Brian W. Montague Report LEP Note 165, CERN, 1979
%
%  See also atx atlinopt
clight = PhysConstant.speed_of_light_in_vacuum.value;   % m/s

[get_chrom, varargs]=getflag(varargin, 'get_chrom');
[get_w, varargs]=getflag(varargs, 'get_w');
[twiss_in,varargs]=getoption(varargs,'twiss_in',[]);
[orbitin,varargs]=getoption(varargs,'orbit',[]);
[DPStep,varargs]=getoption(varargs,'DPStep');
[XYStep,varargs]=getoption(varargs,'XYStep');
[dp,varargs]=getoption(varargs,'dp',NaN);
[refpts,varargs]=getargs(varargs,1);

lgth = findspos(ring,length(ring)+1);
if isempty(twiss_in)        % Circular machine
    is6d=check_radiation(ring);
    orbitin=find_orbit(ring,dp,'orbit',orbitin,varargs{:});
    [vps,ms,orbs,mu,ri,ai]=build_1turn_map(ring,dp,refpts,orbitin);
else                        % Transfer line
    if isempty(orbitin), orbitin=zeros(6,1); end
    sigma=build_sigma(twiss_in);
    nv=size(sigma,1);
    dms=nv/2;
    is6d=(dms>=3);
    [vps,ms,orbs,mu,ri,ai]=build_1turn_map(ring,dp,refpts,orbitin,'mxx',sigma*jmat(dms));
end

tunes=mod(angle(vps)/2/pi,1);
damping_rates=-log(abs(vps));
damping_times=lgth / clight ./ damping_rates;

elemdata=struct('R',ri,'M',ms,'A',ai,...
    'mu',num2cell(unwrap(cat(1,mu{:})),2),...
    'SPos',num2cell(findspos(ring,refpts))',...
    'ClosedOrbit',num2cell(orbs,1)'...
    );
ringdata=struct('tune',tunes,'damping_time',damping_times);

if is6d                     % 6D processing
    [alpha,beta,disp]=cellfun(@output6,ri,'UniformOutput',false);
    [elemdata.alpha]=deal(alpha{:});
    [elemdata.beta]=deal(beta{:});
    [elemdata.Dispersion]=deal(disp{:});
    if get_w || get_chrom
        frf=get_rf_frequency(ring);
        DFStep=-DPStep*mcf(atradoff(ring))*frf;
        rgup=atsetcavity(ring,'Frequency',frf+0.5*DFStep);
        rgdn=atsetcavity(ring,'Frequency',frf-0.5*DFStep);
        if get_w
            [ringdata.chromaticity,~,w]=get_disp(rgup,rgdn,NaN,NaN,refpts);
            [elemdata.W]=deal(w{:});
        else
            [ringdata.chromaticity,~,~]=get_disp(rgup,rgdn,NaN,NaN,[]);
        end
    end
else                        % 4D processing
    dp=orbitin(5);
    [alpha,beta]=cellfun(@output4,ri,'UniformOutput',false);
    [elemdata.alpha]=deal(alpha{:});
    [elemdata.beta]=deal(beta{:});
    if get_w
            [ringdata.chromaticity,disp,w]=get_disp(ring,ring,dp+0.5*DPStep,dp-0.5*DPStep,refpts);
            [elemdata.Dispersion]=deal(disp{:});
            [elemdata.W]=deal(w{:});
    elseif get_chrom
            [ringdata.chromaticity,disp,~]=get_disp(ring,ring,dp+0.5*DPStep,dp-0.5*DPStep,refpts,'skip');
            [elemdata.Dispersion]=deal(disp{:});
    end
end

    function [alpha,beta,disp]=output6(ri)
        % Extract output parameters from R matrices
        alpha= -[ri(2,1,1) ri(4,3,2)];
        beta=[ri(1,1,1) ri(3,3,2)];
        disp=ri(1:4,5,3)/ri(5,5,3);
    end

    function [alpha,beta]=output4(ri)
        % Extract output parameters from R matrices
        alpha= -[ri(2,1,1) ri(4,3,2)];
        beta=[ri(1,1,1) ri(3,3,2)];
    end

    function up = unwrap(p)
        % Unwrap negative jumps in betatron
        jumps = diff([zeros(1,size(p,2));p],1,1) < -1.e-3;
        up = p+cumsum(jumps)*2*pi;
    end

    function orbitin=find_orbit(ring,dp,varargin)
        % Find  the closed orbit in 4D or 6D
        [orbitin,vargs]=getoption(varargin,'orbit',[]);
        [ct,~]=getoption(vargs,'ct',NaN);
        if is6d
            if isfinite(dp),warning('AT:linopt','In 6D, "dp" and "ct" are ignored'); end
            if isempty(orbitin)
                [~,orbitin]=findorbit6(ring,'XYStep',XYStep,'DPStep',DPStep);
            end
        else
            if ~isempty(orbitin)
                if length(orbitin) >= 5
                    dp=orbitin(5);
                end
                orbitin=[orbitin(1:4);dp;0];
            elseif isnan(ct)
                if isnan(dp), dp=0; end
                [~,orbitin]=findorbit4(ring,dp,'XYStep',XYStep);
            else
                [~,orbitin]=findsyncorbit(ring,ct,'XYStep',XYStep);
            end
        end
    end

    function [vals,ms,orbs,phis,rmats,as]=build_1turn_map(ring,dp,refpts,orbit,varargin)
        % Build the initial distribution at entrance of the transfer line
        [skip,vargs]=getflag(varargin,'skip');
        if skip
            if is6d
                [mt,ms]=findm66(ring,[],'orbit',orbit,'XYStep',XYStep,'DPStep',DPStep);
            else
                [mt,ms]=findm44(ring,dp,[],'orbit',orbit,'XYStep',XYStep);
            end
            orbs=linepass(ring,orbit,refpts,'KeepLattice');
            orbs=orbs(1:4,:);
        else
            if is6d
                [mt,ms,orbs]=findm66(ring,refpts,'orbit',orbit,'XYStep',XYStep,'DPStep',DPStep);
            else
                [mt,ms,orbs]=findm44(ring,dp,refpts,'orbit',orbit,'XYStep',XYStep);
            end
        end
        ms=squeeze(num2cell(ms,[1 2]));
        [mxx,~]=getoption(vargs,'mxx',mt);
        [amat0,vals]=amat(mxx);
        [~,phis,rmats,as]=r_analysis(amat0,ms);
    end

    function sigma=build_sigma(twiss_in)
        % build the initial conditions for a transfer line: sigma = R1 + R2
        if isfield(twiss_in,'R')
            sigma=sum(twiss_in.R,3);
        else
            slices=num2cell(reshape(1:4,2,2),1);
            v=num2cell(cat(1,twiss_in.alpha,twiss_in.beta),1);
            sigma=zeros(4,4);
            cellfun(@sigma22,v,slices,'UniformOutput',false);
        end
        
        function sigma22(ab,slc)
            alp=ab(1);
            bet=ab(2);
            sigma(slc,slc)=[bet -alp;-alp (1+alp*alp)/bet];
        end
    end
    
    function f=get_rf_frequency(ring)
        % Get the initial RF frequency
        cavities=ring(atgetcells(ring, 'Frequency'));
        freqs=atgetfieldvalues(cavities,'Frequency');
        f=freqs(1);
    end

    function w=chromfunc(ddp,rup,rdn)
        % Compute the chromatic W function
        [aup,bup]=output4(rup);
        [adn,bdn]=output4(rdn);
        db = (bup - bdn) / ddp;
        mb = (bup + bdn) / 2;
        da = (aup - adn) / ddp;
        ma = (aup + adn) / 2;
        w = sqrt((da - ma ./ mb .* db).^2 + (db ./ mb).^2);
    end

    function [chrom,disp,w]=get_disp(ringup,ringdn,dpup,dpdn,refpts,varargin)
        % Compute chromaticity, dispersion and W
        [dpup,tuneup,orbup,rup]=offmom(ringup,dpup,refpts);
        [dpdn,tunedn,orbdn,rdn]=offmom(ringdn,dpdn,refpts);
        deltap=dpup-dpdn;
        chrom=(tuneup-tunedn)./deltap;
        disp=num2cell((orbup-orbdn)./deltap,1)';
        w=cellfun(@(r1,r2) chromfunc(deltap,r1,r2),rup,rdn,'UniformOutput',false);
        
        function [dp,tunes,orbs,rmats]=offmom(ring,dp,refpts)
            orbit=find_orbit(ring,dp);
            [vals,~,orbs,~,rmats,~]=build_1turn_map(ring,dp,refpts,orbit,varargin{:});
            tunes=mod(angle(vals)/2/pi,1);
            dp=orbit(5);
        end
    end

end

