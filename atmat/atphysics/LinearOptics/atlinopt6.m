function [ringdata,elemdata] = atlinopt6(ring, varargin)
%ATLINOPT6 Performs linear analysis of the lattice
%
% [RINGDATA,ELEMDATA] = ATLINOPT6(RING,REFPTS)
%
%For circular machines, ATLINOPT6 analyses
%the 4x4 1-turn transfer matrix if radiation is OFF, or
%the 6x6 1-turn transfer matrix if radiation is ON.
%
%For a transfer line, The "twiss_in" intput must contain either:
% - a field 'R', as provided by ATLINOPT6, or
% - the fields 'beta' and 'alpha', as provided by ATLINOPT and ATLINOPT6
%
%RINGDATA is a structure array with fields:
%   tune            Fractional tunes
%   damping_time    Damping times [s]
%   chromaticity    Chromaticities
%
%ELEMDATA is a structure array with fields:
%   SPos        - longitudinal position [m]
%   ClosedOrbit - 6x1 closed orbit vector with components
%                 x, px, y, py, dp/d, ct (momentums, NOT angles)
%   Dispersion  - [eta_x; eta'_x; eta_y; eta'_y] 4x1 dispersion vector
%   R           - DxDx(D/2) R-matrices
%   A           - DxD A-matrix
%   M           - DxD transfer matrix M from the beginning of RING
%   beta        - [betax, betay]                 1x2 beta vector
%   alpha       - [alphax, alphay]               1x2 alpha vector
%   mu          - [mux, muy] 	Betatron phase advances
%   W           - [Wx, Wy]  Chromatic amplitude function [3] (only with the
%                           get_w flag)
%
%   Use the Matlab function "cat" to get the data from fields of ELEMDATA as MATLAB arrays.
%   Example:
%   >> [ringdata, elemdata] = ATLINOPT6(ring,1:length(ring));
%   >> beta = cat(1,elemdata.beta);
%   >> s = cat(1,elemdata.SPos);
%   >> plot(S,beta)
%
%   All values are specified at the entrance of each element specified in REFPTS.
%   REFPTS is an array of increasing indexes that  select elements
%   from the range 1 to length(LINE)+1. Defaults to 1 (initial point)
%   See further explanation of REFPTS in the 'help' for FINDSPOS
%
% [...] = ATLINOPT6(...,'get_chrom')
%   Trigger the computation of chromaticities
%
% [...] = ATLINOPT6(...,'get_w')
%   Trigger the computation of chromatic amplitude functions (time consuming)
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
% [...] = ATLINOPT6(...,'dp',DP)
%   Specify off-momentum 
%
% [...] = ATLINOPT6(...,'dct',DCT)
%   Specify path lengthening
%
%  REFERENCES
%   [1] Etienne Forest, Phys. Rev. E 58, 2481 – Published 1 August 1998
%   [2] Andrzej Wolski, Phys. Rev. ST Accel. Beams 9, 024001 –
%       Published 3 February 2006
%   [3] Brian W. Montague Report LEP Note 165, CERN, 1979
%
%  See also atlinopt atlinopt2 atlinopt4 tunechrom

[ringdata,elemdata]=wrapper6d(ring,@xlinopt6,varargin{:});

    function [ringdata,elemdata] = xlinopt6(ring, is6d, varargin)
        clight = PhysConstant.speed_of_light_in_vacuum.value;   % m/s
        [get_chrom, varargs]=getflag(varargin, 'get_chrom');
        [get_w, varargs]=getflag(varargs, 'get_w');
        [dpargs,varargs]=getoption(varargs,{'orbit','dp','dct'});
        [twiss_in,varargs]=getoption(varargs,'twiss_in',[]);
        [DPStep,~]=getoption(varargs,'DPStep');
        [refpts,varargs]=getargs(varargs,1);
        [cavpts,varargs]=getoption(varargs,'cavpts',[]);

        if isempty(twiss_in)        % Circular machine
            [orbs,orbitin]=findorbit(ring,refpts,dpargs{:},varargs{:});
            [vps,ms,mu,ri,ai]=build_1turn_map(ring,refpts,orbitin);
            o1P=[];
            o1M=[];
        else                        % Transfer line
            [orbitin,sigma,d0]=build_sigma(twiss_in,getoption(dpargs,'orbit',[]));
            orbs=linepass(ring,orbitin,refpts);
            dorbit=0.5*DPStep*[d0;1;0];
            nv=size(sigma,1);
            dms=nv/2;
            [vps,ms,mu,ri,ai]=build_1turn_map(ring,refpts,orbitin,'mxx',sigma*jmat(dms));
            o1P=orbitin+dorbit;
            o1M=orbitin-dorbit;
        end

        lgth = findspos(ring,length(ring)+1);
        tunes=mod(angle(vps)/2/pi,1);
        damping_rates=-log(abs(vps));
        damping_times=lgth / clight ./ damping_rates;

        elemdata=struct('R',ri,'M',ms,'A',ai,...
            'mu',num2cell(unwrap(cat(1,mu{:})),2),...
            'SPos',num2cell(findspos(ring,refpts))',...
            'ClosedOrbit',num2cell(orbs,1)'...
            )';
        ringdata=struct('tune',tunes,'damping_time',damping_times);

        if is6d                     % 6D processing
            [alpha,beta,disp]=cellfun(@output6,ri,'UniformOutput',false);
            if get_w || get_chrom
                frf=get_rf_frequency(ring);
                DFStep=-DPStep*mcf(atradoff(ring))*frf;
                rgup=atsetcavity(ring,'Frequency',frf+0.5*DFStep,'cavpts',cavpts);
                rgdn=atsetcavity(ring,'Frequency',frf-0.5*DFStep,'cavpts',cavpts);
                [~,o1P]=findorbit(rgup,'guess',orbitin,varargs{:});
                [~,o1M]=findorbit(rgdn,'guess',orbitin,varargs{:});
                if get_w
                    [ringdata.chromaticity,w]=chrom_w(rgup,rgdn,o1P,o1M,refpts);
                    [elemdata.W]=deal(w{:});
                else
                    tuneP=tune6(rgup,'orbit',o1P,varargs{:});
                    tuneM=tune6(rgdn,'orbit',o1M,varargs{:});
                    deltap=o1P(5)-o1M(5);
                    ringdata.chromaticity = (tuneP - tuneM)/deltap;
                end
            end
        else                        % 4D processing
            dp=orbitin(5);
            [alpha,beta]=cellfun(@output4,ri,'UniformOutput',false);
            [orbitP,o1P]=findorbit4(ring,dp+0.5*DPStep,refpts,'orbit',o1P,'guess',orbitin,varargs{:});
            [orbitM,o1M]=findorbit4(ring,dp-0.5*DPStep,refpts,'orbit',o1M,'guess',orbitin,varargs{:});
            disp = num2cell((orbitP-orbitM)/DPStep,1);
            if get_w
                [ringdata.chromaticity,w]=chrom_w(ring,ring,o1P,o1M,refpts);
                [elemdata.W]=deal(w{:});
            elseif get_chrom
                tuneP=tune4(ring,'orbit',o1P,varargs{:});
                tuneM=tune4(ring,'orbit',o1M,varargs{:});
                ringdata.chromaticity = (tuneP - tuneM)/DPStep;
            end
        end

        [elemdata.alpha]=deal(alpha{:});
        [elemdata.beta]=deal(beta{:});
        [elemdata.Dispersion]=deal(disp{:});

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

        function [vals,ms,phis,rmats,as]=build_1turn_map(ring,refpts,orbit,varargin)
            % Build the initial distribution at entrance of the transfer line
            if is6d
                [mt,ms]=findm66(ring,refpts,'orbit',orbit,varargs{:});
            else
                [mt,ms]=findm44(ring,NaN,refpts,'orbit',orbit,varargs{:});
            end
            ms=squeeze(num2cell(ms,[1 2]));
            [mxx,~]=getoption(varargin,'mxx',mt);
            [amat0,vals]=amat(mxx);
            [phis,rmats,as]=r_analysis(amat0,ms);
        end

        function [orbit,sigma,d0]=build_sigma(twiss_in,orbit)
            % build the initial conditions for a transfer line: sigma = R1 + R2
            if isempty(orbit)
                orbit=zeros(6,1);
                if isfield(twiss_in, 'ClosedOrbit')
                    lo=length(twiss_in.ClosedOrbit);
                    orbit(1:lo)=twiss_in.ClosedOrbit;
                end
            end
            if isfield(twiss_in,'R')
                sz=size(twiss_in.R);
                % For some reason, "emittances" must be different...
                sigma=twiss_in.R(:,:,1)+10.0*twiss_in.R(:,:,2);
                if sz(3) >= 3
                    sigma=sigma+0.1*twiss_in.R(:,:,3);
                end
            else
                slices=num2cell(reshape(1:4,2,2),1);
                v=num2cell(cat(1,twiss_in.alpha,twiss_in.beta),1);
                sigma=zeros(4,4);
                cellfun(@sigma22,v,slices,'UniformOutput',false);
            end
            if isfield(twiss_in,'Dispersion')
                d0=twiss_in.Dispersion;
            else
                d0=zeros(4,1);
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

        function [chrom,w]=chrom_w(ringup,ringdn,orbup,orbdn,refpts)
            % Compute chromaticity, dispersion and W
            [dpup,tuneup,rup]=offmom(ringup,orbup,refpts);
            [dpdn,tunedn,rdn]=offmom(ringdn,orbdn,refpts);
            delp=dpup-dpdn;
            chrom=(tuneup-tunedn)./delp;
            w=cellfun(@(r1,r2) chromfunc(delp,r1,r2),rup,rdn,'UniformOutput',false);

            function [dp,tunes,rmats]=offmom(ring,orbit,refpts)
                dp=orbit(5);
                [vals,~,~,rmats,~]=build_1turn_map(ring,refpts,orbit);
                tunes=mod(angle(vals)/2/pi,1);
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
        end

        function tune=tune4(ring,varargin)
            mm=findm44(ring,NaN,varargin{:});
            [~,vals]=amat(mm);
            tune=mod(angle(vals)/2/pi,1);
        end

        function tune=tune6(ring,varargin)
            mm=findm66(ring,varargin{:});
            [~,vals]=amat(mm);
            tune=mod(angle(vals)/2/pi,1);
        end
    end
end
